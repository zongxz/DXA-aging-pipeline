# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 01_imaging_model
# 5 折交叉验证训练 ResNet101，输出各部位预测年龄与折外预测
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import gc
import json
import math
import time
import random
import argparse

import numpy as np
import pandas as pd
from PIL import Image

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torchvision import models, transforms

from tqdm import tqdm
from sklearn.metrics import (mean_absolute_error, mean_squared_error, r2_score)
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Configuration
PROJECT_ROOT = "/path/to/project"
METADATA_DIR = os.path.join(PROJECT_ROOT, "data", "model_metadata")
RESULT_DIR = os.path.join(PROJECT_ROOT, "results", "training")

# Supported modes
MODES = [
    "single_left_knee",
    "single_right_knee",
    "single_left_hip",
    "single_right_hip",
    "double_knee",
    "double_hip",
    "double_left",
    "double_right",
    "four_joint",
    "single_right_knee_57047",
    "single_right_hip_57047",
    "double_right_51638",
    "single_right_knee_51638",
    "single_right_hip_51638"
]

# Reproducibility
SEED = 42
NUM_FOLDS = 5
MAX_EPOCHS = 100
EARLY_STOPPING_PATIENCE = 12
MIN_DELTA = 0.0001


# Image
IMAGE_SIZE = 224

# Batch size
#
# ResNet101 + multi-view
# A800 80GB
DEFAULT_BATCH_SIZE = {
    "single_left_knee": 256,
    "single_right_knee": 256,
    "single_left_hip": 256,
    "single_right_hip": 256,
    "single_right_knee_57047": 256,
    "single_right_hip_57047": 256,
    "double_knee": 128,
    "double_hip": 128,
    "double_left": 128,
    "double_right": 128,
    "four_joint": 64,
    "double_right_51638": 128,
    "single_right_knee_51638": 256,
    "single_right_hip_51638": 256    
}

# DataLoader
NUM_WORKERS = 16
PIN_MEMORY = True
PERSISTENT_WORKERS = True


# Optimization
LEARNING_RATE = 1e-4
WEIGHT_DECAY = 1e-4
MIN_LR = 1e-6
LR_PATIENCE = 4
LR_FACTOR = 0.5

# Model
DROPOUT = 0.3

# Seed
def seed_everything(seed=SEED):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    # Reproducibility
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.deterministic = False


# Device
def get_device():
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print("Using GPU:",  torch.cuda.get_device_name(0))
    else:
        device = torch.device("cpu")
        print("Using CPU")
    return device

# Metrics
def calculate_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=np.float64)
    y_pred = np.asarray(y_pred, dtype=np.float64)

    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)

    if (len(y_true) > 1 and np.std(y_true) > 0 and np.std(y_pred) > 0):
        pearson_r, pearson_p = pearsonr(y_true, y_pred)
    else:
        pearson_r = np.nan
        pearson_p = np.nan
        
    return {
        "MAE": float(mae),
        "RMSE": float(rmse),
        "R2": float(r2),
        "Pearson_r": float(pearson_r),
        "Pearson_p": float(pearson_p)
    }

# Image transforms
def get_train_transform():
    return transforms.Compose([
        transforms.Resize((IMAGE_SIZE, IMAGE_SIZE)),
        transforms.RandomHorizontalFlip(p=0.5),
        transforms.RandomRotation(degrees=5),
        transforms.ColorJitter(brightness=0.15, contrast=0.15),
        transforms.ToTensor(),
        transforms.Normalize(
            mean=[0.485, 0.456, 0.406],
            std=[0.229, 0.224, 0.225]
        )
    ])


def get_val_transform():
    return transforms.Compose([
        transforms.Resize((IMAGE_SIZE, IMAGE_SIZE)),
        transforms.ToTensor(),
        transforms.Normalize(
            mean = [0.485, 0.456, 0.406],
            std = [0.229, 0.224, 0.225]
        )
    ])

# Dataset
class XRayAgeDataset(Dataset):
    def __init__(self, df, image_cols, transform=None):
        self.df = df.reset_index(drop=True)
        self.image_cols = image_cols
        self.transform = transform

        # Validate metadata BEFORE training
        for col in self.image_cols:
            if col not in self.df.columns:
                raise ValueError(f"Missing column: {col}")

        # Validate all image paths
        for col in self.image_cols:
            invalid = ~self.df[col].apply(lambda x: isinstance(x, str) and os.path.isfile(x))

            if invalid.any():
                bad_rows = self.df.loc[invalid, ["eid", col]].head()
                
                raise FileNotFoundError(
                    f"Invalid image path "
                    f"in column {col}\n"
                    f"{bad_rows}")

    def __len__(self):
        return len(self.df)

    def load_image(self, path):
        try:
            image = Image.open(path)
            
            # X-ray is grayscale,
            # convert to RGB for ImageNet
            image = image.convert("RGB")
            if self.transform is not None:
                image = self.transform(image)
            return image

        except Exception as e:
            raise RuntimeError(
                f"Cannot read image:\n"
                f"{path}\n"
                f"Error: {e}")


    def __getitem__(self, index):

        row = self.df.iloc[index]
        images = []
        for col in self.image_cols:
            image = self.load_image(row[col])
            images.append(image)

        # [N_views, C, H, W]
        images = torch.stack(images, dim=0)
        age = torch.tensor(float(row["age"]), dtype=torch.float32)

        return {
            "images": images,
            "age": age,
            "eid": int(row["eid"])
        }

# ResNet101 backbone
def build_resnet101():
    try:
        weights = (models.ResNet101_Weights.IMAGENET1K_V2)
        model = models.resnet101(weights=weights)
    except Exception:
        print("Could not load IMAGENET1K_V2.")
        print("Trying default pretrained weights...")
        model = models.resnet101(pretrained=True)

    feature_dim = (model.fc.in_features)
    model.fc = nn.Identity()
    
    return model, feature_dim

# Multi-view ResNet101
class MultiViewResNet101(nn.Module):

    def __init__(self, num_views, dropout=DROPOUT):
        super().__init__()
        self.num_views = num_views
        
        # One shared ResNet101
        self.backbone, feature_dim = (build_resnet101())

        # Regression head
        total_dim = (feature_dim * num_views)

        self.regressor = nn.Sequential(
            nn.Linear(total_dim, 1024),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(1024, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(256, 1)
        )


    def forward(self, x):
        # x:
        # [B, V, C, H, W]
        batch_size = x.shape[0]
        features = []
        
        # Process each view separately
        #
        # This reduces peak GPU memory compared with
        # flattening B*V directly.
        for view_idx in range(self.num_views):
            view = x[:, view_idx, :, :, :]
            feature = self.backbone(view)
            features.append(feature)
            
        # [B, V*2048]
        features = torch.cat(features, dim=1)
        output = self.regressor(features)
        return output.squeeze(1)

# Train one epoch
def train_one_epoch(model, loader, optimizer, criterion, scaler, device, epoch, total_epochs):

    model.train()
    running_loss = 0.0
    n_samples = 0
    progress = tqdm(loader, desc=(f"Train {epoch}/{total_epochs}"), leave=False)

    for batch in progress:
        images = batch["images"].to(device, non_blocking=True)
        age = batch["age"].to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
    
        # AMP
        with torch.cuda.amp.autocast(enabled=device.type == "cuda"):
            pred = model(images)
            loss = criterion(pred, age)
            
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        
        batch_size = age.size(0)
        running_loss += (loss.item() * batch_size)
        n_samples += (batch_size)

        progress.set_postfix(loss=f"{loss.item():.4f}")

    epoch_loss = (running_loss / max(n_samples, 1))
    
    return epoch_loss


# Validation
@torch.no_grad()
def validate(model, loader, criterion, device):
    model.eval()
    running_loss = 0.0
    n_samples = 0
    all_true = []
    all_pred = []
    all_eid = []

    for batch in tqdm(loader, desc="Validation", leave=False):

        images = batch["images"].to(device, non_blocking=True)
        age = batch["age"].to(device, non_blocking=True)
        with torch.cuda.amp.autocast(enabled=device.type == "cuda"):
            pred = model(images)
            loss = criterion(pred, age)

        batch_size = age.size(0)
        running_loss += (loss.item() * batch_size)
        n_samples += (batch_size)

        all_true.extend(age.cpu().numpy())
        all_pred.extend(pred.cpu().numpy())
        all_eid.extend(batch["eid"].cpu().numpy())

    val_loss = (running_loss / max(n_samples, 1))
    metrics = calculate_metrics(all_true, all_pred)

    return (val_loss, metrics, np.asarray(all_eid), np.asarray(all_true), np.asarray(all_pred))


# Save training curves
def save_training_curves(log_df, output_dir, fold):
    # Loss
    plt.figure(figsize=(8, 6))
    
    plt.plot(log_df["epoch"], log_df["train_loss"], label="Train loss")
    plt.plot(log_df["epoch"], log_df["val_loss"], label="Validation loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title(f"Fold {fold} Loss")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "loss_curve.png"), dpi=300)
    plt.close()

    # MAE
    plt.figure(figsize=(8, 6))
    plt.plot(log_df["epoch"],log_df["train_mae"], label="Train MAE")
    plt.plot(log_df["epoch"], log_df["val_mae"], label="Validation MAE")
    plt.xlabel("Epoch")
    plt.ylabel("MAE (years)")
    plt.title( f"Fold {fold} MAE")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mae_curve.png"), dpi=300)
    plt.close()

# Save OOF prediction
def save_oof(eid, age, pred, sex, fold, output_path):
    df = pd.DataFrame({
        "eid":  eid.astype(np.int64),
        "age": age.astype(float),
        "sex": sex.astype(int),
        "fold": fold.astype(int),
        "predicted_age": pred.astype(float)
    })
    df["age_gap"] = (df["predicted_age"] - df["age"])
    
    # Ensure exactly one prediction per eid
    if df["eid"].duplicated().any():
        raise RuntimeError("OOF contains duplicate eid.")

    df.to_csv(output_path, index=False)
    return df

# Train one fold
def train_one_fold(mode, fold, full_df, image_cols, device, args):
    print("\n")
    print("=" * 100)
    print(f"{mode} | Fold {fold}")
    print("=" * 100)
    
    # Split using PREDEFINED metadata fold
    train_df = full_df[full_df["fold"] != fold].copy()
    val_df = full_df[full_df["fold"] == fold].copy()
    
    print("Train participants:", len(train_df))
    print("Validation participants:", len(val_df))
    
    # Dataset
    train_dataset = XRayAgeDataset(train_df, image_cols, transform=get_train_transform())
    val_dataset = XRayAgeDataset(val_df, image_cols, transform=get_val_transform())

    batch_size = (
        args.batch_size
        if args.batch_size is not None
        else DEFAULT_BATCH_SIZE[mode]
    )

    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=args.num_workers,
        pin_memory=PIN_MEMORY,
        persistent_workers=(PERSISTENT_WORKERS if args.num_workers > 0 else False), 
        drop_last=True
    )

    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=PIN_MEMORY,
        persistent_workers=(PERSISTENT_WORKERS if args.num_workers > 0 else False),
        drop_last=False
    )

    # Model
    num_views = len(image_cols)
    model = MultiViewResNet101(num_views=num_views)
    model = model.to(device)
    
    # Loss
    criterion = nn.SmoothL1Loss(beta=1.0)

    # Optimizer
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)

    # LR scheduler
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=LR_FACTOR, patience=LR_PATIENCE, min_lr=MIN_LR)

    # AMP
    scaler = torch.cuda.amp.GradScaler(enabled=(device.type == "cuda"))

    # Output directory
    fold_dir = os.path.join(RESULT_DIR, mode, f"fold_{fold}")
    os.makedirs(fold_dir, exist_ok=True)

    # Training log
    log = []
    best_val_mae = np.inf
    best_epoch = 0
    patience_counter = 0
    start_time = time.time()
    
    # Training loop
    for epoch in range(1, MAX_EPOCHS + 1):
        current_lr = (optimizer.param_groups[0]["lr"])
        print()
        print(f"Fold {fold} | Epoch {epoch}/{MAX_EPOCHS}")
        print(f"Learning rate: {current_lr:.8f}")

        # Train
        train_loss = train_one_epoch(model, train_loader, optimizer, criterion, scaler, device, epoch, MAX_EPOCHS)

        # Validation
        (val_loss, val_metrics, val_eid, val_age, val_pred) = validate(model, val_loader, criterion, device)

        # Train metrics
        #
        # To avoid another complete training pass,
        # use train loss here.

        # SmoothL1 loss is the optimization loss.
        # For train MAE, we calculate it approximately
        # from the validation-independent training loader
        # only if requested.
        #
        # For speed, train_mae is recorded as NaN.
        # Validation MAE is the key metric for checkpointing.
        train_mae = np.nan

        # Scheduler
        scheduler.step(val_metrics["MAE"])

        # Log
        row = {
            "epoch": epoch,
            "learning_rate": current_lr,
            "train_loss": train_loss,
            "train_mae": train_mae,
            "val_loss": val_loss,
            "val_mae": val_metrics["MAE"],
            "val_rmse": val_metrics["RMSE"],
            "val_r2": val_metrics["R2"],
            "val_pearson_r": val_metrics["Pearson_r"],
            "val_pearson_p": val_metrics["Pearson_p"]
        }
        log.append(row)
        print()
        print(f"Train loss: {train_loss:.5f}")
        print(f"Val loss: {val_loss:.5f}")
        print(f"Val MAE: {val_metrics['MAE']:.5f}")
        print(f"Val RMSE: {val_metrics['RMSE']:.5f}")
        print(f"Val R²: {val_metrics['R2']:.5f}")
        print(f"Val Pearson r: {val_metrics['Pearson_r']:.5f}")

        # Save log every epoch
        log_df = pd.DataFrame(log)
        log_df.to_csv(os.path.join(fold_dir, "training_log.csv"), index=False)

        # Best model
        if (val_metrics["MAE"] < best_val_mae - MIN_DELTA):
            best_val_mae = (val_metrics["MAE"])
            best_epoch = epoch
            patience_counter = 0
            checkpoint = {
                "epoch": epoch,
                "mode": mode,
                "fold": fold,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "scheduler_state_dict": scheduler.state_dict(),
                "best_val_mae": best_val_mae,
                "image_cols": image_cols,
                "image_size": IMAGE_SIZE,
                "num_views": num_views
            }

            torch.save(checkpoint, os.path.join(fold_dir, "best_model.pth"))
            print()
            print("NEW BEST MODEL")
            print(f"Best Val MAE: {best_val_mae:.5f}")

        else:
            patience_counter += 1

        # Curves
        save_training_curves(log_df, fold_dir, fold)

        # Early stopping
        if (patience_counter >= EARLY_STOPPING_PATIENCE):
            print()
            print("Early stopping.")
            print(f"Best epoch: {best_epoch}")
            print(f"Best Val MAE: {best_val_mae:.5f}")
            break

    # Load BEST checkpoint
    checkpoint_path = os.path.join(fold_dir, "best_model.pth")
    checkpoint = torch.load(checkpoint_path, map_location=device)

    model.load_state_dict(checkpoint["model_state_dict"])

    # Final OOF prediction
    # using BEST model
    (final_val_loss, final_metrics, final_eid, final_age, final_pred) = validate(model, val_loader, criterion, device)

    # Sex
    sex_map = (val_df.set_index("eid")["sex"].to_dict())
    final_sex = np.asarray([sex_map[int(eid)] for eid in final_eid])

    # Save fold prediction
    fold_prediction = pd.DataFrame({
        "eid": final_eid.astype(np.int64),
        "age": final_age.astype(float),
        "sex": final_sex.astype(int),
        "fold": np.full(len(final_eid), fold, dtype=int),
        "predicted_age": final_pred.astype(float)
    })
    fold_prediction["age_gap"] = (fold_prediction["predicted_age"] - fold_prediction["age"])
    fold_prediction.to_csv(os.path.join(fold_dir, "fold_oof.csv"),index=False)

    # Fold metrics
    fold_metrics = {
        "mode": mode,
        "fold": fold,
        "N": len(final_eid),
        "best_epoch": best_epoch,
        "best_val_mae": best_val_mae,
        "MAE": final_metrics["MAE"],
        "RMSE": final_metrics["RMSE"],
        "R2": final_metrics["R2"],
        "Pearson_r": final_metrics["Pearson_r"],
        "Pearson_p": final_metrics["Pearson_p"]
    }

    with open(os.path.join(fold_dir, "fold_metrics.json"), "w") as f:
        json.dump(fold_metrics, f, indent=4)

    print()
    print(f"Fold {fold} finished.")
    print(f"Best epoch: {best_epoch}")
    print(f"OOF MAE: {final_metrics['MAE']:.5f}")
    print(f"OOF RMSE: {final_metrics['RMSE']:.5f}")

    # Cleanup
    del model
    del optimizer
    del scheduler
    del scaler
    del train_loader
    del val_loader
    del train_dataset
    del val_dataset

    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    elapsed = (time.time() - start_time)
    print(f"Fold time: {elapsed / 3600:.2f} hours")

    return (fold_prediction, fold_metrics)

# Train one mode
def train_one_mode(mode, device, args):
    print("\n")
    print("#" * 100)
    print(f"MODE: {mode}")
    print("#" * 100)

    # Load metadata
    metadata_path = os.path.join(METADATA_DIR,  mode + ".csv")
    
    if not os.path.isfile(metadata_path):
        raise FileNotFoundError(f"Metadata not found:\n{metadata_path}")
        
    df = pd.read_csv(metadata_path)

    # Validate metadata
    required = ["eid", "age", "sex", "fold"]

    for col in required:
        if col not in df.columns:
            raise ValueError(f"{mode}: missing {col}")

    image_cols = [c for c in df.columns if c.startswith("image")]

    if len(image_cols) == 0:
        raise ValueError(f"{mode}: No image columns.")

    # Basic checks
    if df["eid"].duplicated().any():
        raise RuntimeError(f"{mode}: duplicate eid.")

    if not set(df["fold"].unique()).issubset(set(range(1, NUM_FOLDS + 1))):
        raise RuntimeError(f"{mode}: Invalid fold values.")

    for col in image_cols:
        invalid = ~df[col].apply(lambda x: isinstance(x, str) and os.path.isfile(x))
        
        if invalid.any():
            raise FileNotFoundError(f"{mode}: invalid images in {col}")
    
    print("Participants:", len(df))
    print("Number of views:", len(image_cols))
    print("Image columns:", image_cols)
    print("Fold counts:")

    print(df["fold"].value_counts().sort_index().to_dict())

    # Store all fold predictions
    all_predictions = []
    all_fold_metrics = []

    # Five folds
    for fold in range(1, NUM_FOLDS + 1):
        fold_prediction, fold_metrics = (
            train_one_fold(mode, fold, df, image_cols, device, args))

        all_predictions.append(fold_prediction)
        all_fold_metrics.append(fold_metrics)

    # Combine OOF
    oof = pd.concat(all_predictions, axis=0, ignore_index=True)

    # Sort by eid
    oof = oof.sort_values("eid").reset_index(drop=True)

    # Critical OOF checks
    if oof["eid"].duplicated().any():
        raise RuntimeError(f"{mode}: OOF has duplicate eid.")

    if len(oof) != len(df):
        raise RuntimeError(f"{mode}: OOF N={len(oof)} but metadata N={len(df)}")

    if set(oof["eid"]) != set(df["eid"]):
        raise RuntimeError(f"{mode}: OOF eid does not exactly match metadata eid.")

    # Overall OOF metrics
    overall_metrics = calculate_metrics(oof["age"].values, oof["predicted_age"].values)

    overall_result = {
        "mode": mode,
        "N": len(oof),
        "MAE": overall_metrics["MAE"],
        "RMSE": overall_metrics["RMSE"],
        "R2": overall_metrics["R2"],
        "Pearson_r": overall_metrics["Pearson_r"],
        "Pearson_p": overall_metrics["Pearson_p"]
    }

    # Save OOF
    mode_dir = os.path.join(RESULT_DIR, mode)
    os.makedirs(mode_dir, exist_ok=True)
    oof_path = os.path.join(mode_dir, f"{mode}_oof.csv")
    oof.to_csv(oof_path,index=False)

    # Save fold metrics
    fold_metrics_df = pd.DataFrame(all_fold_metrics)
    fold_metrics_df.to_csv(os.path.join(mode_dir,f"{mode}_fold_metrics.csv"),index=False)

    # Save overall metrics
    overall_df = pd.DataFrame([overall_result])
    overall_df.to_csv(os.path.join(mode_dir,f"{mode}_overall_metrics.csv"),index=False)

    # Print
    print("\n")
    print("=" * 100)
    print(f"{mode} FINAL OOF RESULT")
    print("=" * 100)
    print(f"N          :{len(oof)}")
    print(f"MAE        :{overall_metrics['MAE']:.5f}")
    print(f"RMSE       :{overall_metrics['RMSE']:.5f}")
    print(f"R2         :{overall_metrics['R2']:.5f}" )
    print(f"Pearson r  :{overall_metrics['Pearson_r']:.5f}")
    print()
    print("OOF:")
    print(oof_path)

    return overall_result

# Main
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str, required=True, choices=MODES)
    parser.add_argument("--batch_size", type=int, default=None)
    parser.add_argument("--num_workers", type=int, default=NUM_WORKERS)
    args = parser.parse_args()

    # Seed
    seed_everything(SEED)
    
    # Output directory
    os.makedirs(RESULT_DIR, exist_ok=True)

    # Device
    device = get_device()

    # GPU information
    if device.type == "cuda":
        print("GPU memory:", round(torch.cuda.get_device_properties(0).total_memory / 1024**3, 2), "GB")

    print()
    print("Mode:", args.mode)
    print("Model: ResNet101")
    print("Image size:", IMAGE_SIZE)
    print("Max epochs:", MAX_EPOCHS)
    print("Early stopping patience:", EARLY_STOPPING_PATIENCE)
    print("Learning rate:", LEARNING_RATE)

    # Train
    result = train_one_mode(args.mode, device, args)

    # Final
    print("\n")
    print("=" * 100)
    print("TRAINING COMPLETE")
    print("=" * 100)
    print(json.dumps(result, indent=4))

# Run
if __name__ == "__main__":
    main()