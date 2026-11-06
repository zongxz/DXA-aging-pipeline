# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 01_imaging_model
# 外部数据集验证：左右部位交叉验证与年龄预测精度评估
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import gc
import json
import argparse

import numpy as np
import pandas as pd

from PIL import Image

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torchvision import models, transforms

from tqdm import tqdm

from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score
)

from scipy.stats import pearsonr

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# Project
# ============================================================

PROJECT_ROOT = (
    "/path/to/project"
)


# ============================================================
# Training results
# ============================================================

TRAINING_DIR = os.path.join(PROJECT_ROOT, "results", "training")


# ============================================================
# External validation output
# ============================================================

EXTERNAL_DIR = os.path.join(PROJECT_ROOT, "results", "external_validation")


# ============================================================
# Image
# ============================================================

IMAGE_SIZE = 224


# ============================================================
# Batch size
# ============================================================

BATCH_SIZE = {
    "single_right_knee_51638": 128,
    "single_right_hip_51638": 128,
    "double_right_51638": 64

}


# ============================================================
# DataLoader
# ============================================================

DEFAULT_NUM_WORKERS = 16


# ============================================================
# Five trained folds
# ============================================================

FOLDS = [1, 2, 3, 4, 5]


# ============================================================
# Model configuration
# ============================================================

MODEL_MODES = {

    "single_right_knee_51638": {
        "num_views": 1,
        "image_cols": ["image1"]
    },

    "single_right_hip_51638": {
        "num_views": 1,
        "image_cols": ["image1"]
    },
    
    "double_right_51638": {
        "num_views": 2,
        "image_cols": ["image1", "image2"]
    }
}


# ============================================================
# Validation transform
#
# MUST be the same as training validation transform
# ============================================================

def get_transform():

    return transforms.Compose([

        transforms.Resize((IMAGE_SIZE, IMAGE_SIZE)),

        transforms.ToTensor(),

        transforms.Normalize(

            mean=[
                0.485,
                0.456,
                0.406
            ],
            std=[0.229, 0.224, 0.225]
        )

    ])


# ============================================================
# Dataset
# ============================================================

class ExternalXRayDataset(
    Dataset
):

    def __init__(
        self,
        df,
        image_cols,
        transform
    ):

        self.df = df.reset_index(
            drop=True
        )

        self.image_cols = image_cols

        self.transform = transform


        # ----------------------------------------------------
        # Check columns
        # ----------------------------------------------------

        required = [
            "eid",
            "age",
            "sex"
        ]

        for col in required:

            if col not in self.df.columns:

                raise ValueError(
                    f"Missing required column: {col}"
                )


        for col in self.image_cols:

            if col not in self.df.columns:

                raise ValueError(
                    f"Missing image column: {col}"
                )


        # ----------------------------------------------------
        # Check duplicate eid
        # ----------------------------------------------------

        if self.df["eid"].duplicated().any():

            duplicated = self.df.loc[
                self.df["eid"].duplicated(),
                "eid"
            ].tolist()[:20]


            raise RuntimeError(

                "Duplicate eid detected.\n"
                f"Examples: {duplicated}"

            )


        # ----------------------------------------------------
        # Check image files
        # ----------------------------------------------------

        for col in self.image_cols:

            invalid = ~self.df[col].apply(

                lambda x:
                isinstance(x, str)
                and
                os.path.isfile(x)

            )


            if invalid.any():

                bad = self.df.loc[
                    invalid,
                    ["eid", col]
                ].head(20)


                raise FileNotFoundError(

                    f"Missing image files "
                    f"in {col}:\n"
                    f"{bad}"

                )


    def __len__(self):

        return len(
            self.df
        )


    def load_image(
        self,
        path
    ):

        try:

            image = Image.open(
                path
            ).convert(
                "RGB"
            )

            image = self.transform(
                image
            )

            return image

        except Exception as e:

            raise RuntimeError(

                f"Failed to load image:\n"
                f"{path}\n"
                f"{e}"

            )


    def __getitem__(
        self,
        index
    ):

        row = self.df.iloc[
            index
        ]


        images = []


        for col in self.image_cols:

            image = self.load_image(
                row[col]
            )

            images.append(
                image
            )


        images = torch.stack(
            images,
            dim=0
        )


        return {

            "images":
                images,

            "age":
                torch.tensor(
                    float(row["age"]),
                    dtype=torch.float32
                ),

            "sex":
                torch.tensor(
                    int(row["sex"]),
                    dtype=torch.long
                ),

            "eid":
                int(row["eid"])

        }


# ============================================================
# ResNet101 model
#
# MUST match training architecture
# ============================================================

class MultiViewResNet101(
    nn.Module
):

    def __init__(
        self,
        num_views,
        dropout=0.3
    ):

        super().__init__()


        self.num_views = num_views


        # ----------------------------------------------------
        # ResNet101
        # ----------------------------------------------------

        self.backbone = models.resnet101(
            weights=None
        )


        feature_dim = (
            self.backbone.fc.in_features
        )


        self.backbone.fc = nn.Identity()


        # ----------------------------------------------------
        # Regression head
        # ----------------------------------------------------

        total_dim = (
            feature_dim
            *
            num_views
        )


        self.regressor = nn.Sequential(

            nn.Linear(
                total_dim,
                1024
            ),

            nn.ReLU(
                inplace=True
            ),

            nn.Dropout(
                dropout
            ),

            nn.Linear(
                1024,
                256
            ),

            nn.ReLU(
                inplace=True
            ),

            nn.Dropout(
                dropout
            ),

            nn.Linear(
                256,
                1
            )

        )


    def forward(
        self,
        x
    ):

        features = []


        for view_idx in range(
            self.num_views
        ):

            view = x[
                :,
                view_idx,
                :, :, :
            ]


            feature = self.backbone(
                view
            )


            features.append(
                feature
            )


        features = torch.cat(
            features,
            dim=1
        )


        output = self.regressor(
            features
        )


        return output.squeeze(
            1
        )


# ============================================================
# Load checkpoint
# ============================================================

def load_model(
    checkpoint_path,
    num_views,
    device
):

    print()

    print(
        "Loading:"
    )

    print(
        checkpoint_path
    )


    if not os.path.isfile(
        checkpoint_path
    ):

        raise FileNotFoundError(
            checkpoint_path
        )


    model = MultiViewResNet101(
        num_views=num_views
    )


    checkpoint = torch.load(

        checkpoint_path,

        map_location=device

    )


    # --------------------------------------------------------
    # Check checkpoint
    # --------------------------------------------------------

    if (
        "model_state_dict"
        not in checkpoint
    ):

        raise RuntimeError(

            f"Invalid checkpoint:\n"
            f"{checkpoint_path}"

        )


    model.load_state_dict(

        checkpoint[
            "model_state_dict"
        ]

    )


    model = model.to(
        device
    )


    model.eval()


    if "epoch" in checkpoint:

        print(
            "Checkpoint epoch:",
            checkpoint["epoch"]
        )


    if "best_val_mae" in checkpoint:

        print(
            "Best validation MAE:",
            checkpoint[
                "best_val_mae"
            ]
        )


    return model


# ============================================================
# Prediction
# ============================================================

@torch.no_grad()
def predict(
    model,
    loader,
    device
):

    model.eval()


    all_eid = []

    all_age = []

    all_sex = []

    all_pred = []


    for batch in tqdm(

        loader,

        desc="Predicting",

        leave=False

    ):

        images = batch[
            "images"
        ].to(

            device,

            non_blocking=True

        )


        with torch.cuda.amp.autocast(

            enabled=(
                device.type == "cuda"
            )

        ):

            pred = model(
                images
            )


        all_eid.extend(

            batch[
                "eid"
            ].numpy()

        )


        all_age.extend(

            batch[
                "age"
            ].numpy()

        )


        all_sex.extend(

            batch[
                "sex"
            ].numpy()

        )


        all_pred.extend(

            pred.cpu().numpy()

        )


    result = pd.DataFrame({

        "eid":
            np.asarray(
                all_eid,
                dtype=np.int64
            ),

        "age":
            np.asarray(
                all_age,
                dtype=float
            ),

        "sex":
            np.asarray(
                all_sex,
                dtype=int
            ),

        "predicted_age":
            np.asarray(
                all_pred,
                dtype=float
            )

    })


    return result


# ============================================================
# Metrics
# ============================================================

def calculate_metrics(
    y_true,
    y_pred
):

    y_true = np.asarray(
        y_true,
        dtype=float
    )

    y_pred = np.asarray(
        y_pred,
        dtype=float
    )


    mae = mean_absolute_error(
        y_true,
        y_pred
    )


    rmse = np.sqrt(

        mean_squared_error(
            y_true,
            y_pred
        )

    )


    r2 = r2_score(
        y_true,
        y_pred
    )


    if (

        len(y_true) > 1

        and

        np.std(y_true) > 0

        and

        np.std(y_pred) > 0

    ):

        r, p = pearsonr(
            y_true,
            y_pred
        )

    else:

        r = np.nan

        p = np.nan


    return {

        "N":
            len(y_true),

        "MAE":
            float(mae),

        "RMSE":
            float(rmse),

        "R2":
            float(r2),

        "Pearson_r":
            float(r),

        "Pearson_p":
            float(p)

    }


# ============================================================
# Scatter plot
# ============================================================

def save_scatter_plot(
    df,
    output_path,
    title
):

    plt.figure(
        figsize=(7, 7)
    )


    plt.scatter(

        df["age"],

        df["predicted_age"],

        s=8,

        alpha=0.4

    )


    min_age = min(

        df["age"].min(),

        df["predicted_age"].min()

    )


    max_age = max(

        df["age"].max(),

        df["predicted_age"].max()

    )


    plt.plot(

        [min_age, max_age],

        [min_age, max_age],

        linestyle="--"

    )


    plt.xlabel(
        "Chronological age"
    )

    plt.ylabel(
        "Predicted age"
    )


    plt.title(
        title
    )


    plt.tight_layout()


    plt.savefig(

        output_path,

        dpi=300

    )


    plt.close()


# ============================================================
# Bland-Altman
# ============================================================

def save_bland_altman(
    df,
    output_path,
    title
):

    age = df["age"].values

    pred = df["predicted_age"].values


    mean_age = (
        age + pred
    ) / 2


    difference = (
        pred - age
    )


    mean_difference = (
        np.mean(difference)
    )


    sd_difference = (
        np.std(
            difference,
            ddof=1
        )
    )


    upper = (
        mean_difference
        +
        1.96
        *
        sd_difference
    )


    lower = (
        mean_difference
        -
        1.96
        *
        sd_difference
    )


    plt.figure(
        figsize=(8, 6)
    )


    plt.scatter(

        mean_age,

        difference,

        s=8,

        alpha=0.4

    )


    plt.axhline(
        mean_difference,
        linestyle="-"
    )


    plt.axhline(
        upper,
        linestyle="--"
    )


    plt.axhline(
        lower,
        linestyle="--"
    )


    plt.xlabel(
        "Mean of chronological and predicted age"
    )


    plt.ylabel(
        "Predicted age - chronological age"
    )


    plt.title(
        title
    )


    plt.tight_layout()


    plt.savefig(

        output_path,

        dpi=300

    )


    plt.close()


# ============================================================
# Age gap plot
# ============================================================

def save_age_gap_plot(
    df,
    output_path,
    title
):

    plt.figure(
        figsize=(8, 6)
    )


    plt.scatter(

        df["age"],

        df["age_gap"],

        s=8,

        alpha=0.4

    )


    plt.axhline(
        0,
        linestyle="--"
    )


    plt.xlabel(
        "Chronological age"
    )


    plt.ylabel(
        "Age gap (predicted - chronological)"
    )


    plt.title(
        title
    )


    plt.tight_layout()


    plt.savefig(

        output_path,

        dpi=300

    )


    plt.close()


# ============================================================
# One external validation model
# ============================================================

def run_external_validation(
    model_mode,
    metadata_path,
    device,
    num_workers
):

    config = MODEL_MODES[
        model_mode
    ]


    num_views = config[
        "num_views"
    ]


    image_cols = config[
        "image_cols"
    ]


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    output_dir = os.path.join(

        EXTERNAL_DIR,

        model_mode

    )


    os.makedirs(

        output_dir,

        exist_ok=True

    )


    # --------------------------------------------------------
    # Read external metadata
    # --------------------------------------------------------

    print()

    print(
        "=" * 90
    )

    print(
        f"External validation: "
        f"{model_mode}"
    )

    print(
        "=" * 90
    )


    print(
        "Metadata:"
    )

    print(
        metadata_path
    )


    df = pd.read_csv(
        metadata_path
    )


    print(
        "N:",
        len(df)
    )


    print(
        "Unique eid:",
        df["eid"].nunique()
    )


    if len(df) != df["eid"].nunique():

        raise RuntimeError(
            "External metadata "
            "contains duplicate eid."
        )


    # --------------------------------------------------------
    # Dataset
    # --------------------------------------------------------

    dataset = ExternalXRayDataset(

        df,

        image_cols,

        transform=get_transform()

    )


    batch_size = BATCH_SIZE[
        model_mode
    ]


    loader = DataLoader(

        dataset,

        batch_size=batch_size,

        shuffle=False,

        num_workers=num_workers,

        pin_memory=True,

        persistent_workers=(
            num_workers > 0
        ),

        drop_last=False

    )


    # --------------------------------------------------------
    # Five fold predictions
    # --------------------------------------------------------

    fold_predictions = []


    for fold in FOLDS:

        checkpoint_path = os.path.join(

            TRAINING_DIR,

            model_mode,

            f"fold_{fold}",

            "best_model.pth"

        )


        model = load_model(

            checkpoint_path,

            num_views,

            device

        )


        pred_df = predict(

            model,

            loader,

            device

        )


        # ----------------------------------------------------
        # Validate eid order
        # ----------------------------------------------------

        if not np.array_equal(

            pred_df["eid"].values,

            df["eid"].values

        ):

            raise RuntimeError(

                f"Fold {fold}: "
                "eid order mismatch."

            )


        pred_col = (
            f"pred_fold{fold}"
        )


        fold_result = pred_df[
            [
                "eid",
                "age",
                "sex",
                "predicted_age"
            ]
        ].copy()


        fold_result = (
            fold_result
            .rename(
                columns={
                    "predicted_age":
                    pred_col
                }
            )
        )


        fold_predictions.append(
            fold_result
        )


        # ----------------------------------------------------
        # Save individual fold prediction
        # ----------------------------------------------------

        fold_result.to_csv(

            os.path.join(

                output_dir,

                f"fold_{fold}_predictions.csv"

            ),

            index=False

        )


        # ----------------------------------------------------
        # Fold metrics
        # ----------------------------------------------------

        fold_metrics = calculate_metrics(

            pred_df["age"].values,

            pred_df["predicted_age"].values

        )


        print()

        print(
            f"Fold {fold}"
        )

        print(
            f"MAE: "
            f"{fold_metrics['MAE']:.5f}"
        )

        print(
            f"RMSE: "
            f"{fold_metrics['RMSE']:.5f}"
        )

        print(
            f"R2: "
            f"{fold_metrics['R2']:.5f}"
        )

        print(
            f"Pearson r: "
            f"{fold_metrics['Pearson_r']:.5f}"
        )


        # ----------------------------------------------------
        # Cleanup
        # ----------------------------------------------------

        del model

        gc.collect()

        if torch.cuda.is_available():

            torch.cuda.empty_cache()

    # Combine five fold predictions
    # --------------------------------------------------------
    # Safer explicit merge
    # --------------------------------------------------------

    ensemble = fold_predictions[0].copy()

    for i in range(
        1,
        len(fold_predictions)
    ):

        ensemble = ensemble.merge(

            fold_predictions[i][
                [
                    "eid",
                    f"pred_fold{i + 1}"
                ]
            ],

            on="eid",

            how="inner"

        )


    # --------------------------------------------------------
    # Check
    # --------------------------------------------------------

    if len(ensemble) != len(df):

        raise RuntimeError(

            "Ensemble N mismatch."

        )


    pred_columns = [

        f"pred_fold{fold}"

        for fold in FOLDS

    ]


    # --------------------------------------------------------
    # 5-fold ensemble
    # --------------------------------------------------------

    ensemble[
        "predicted_age"
    ] = ensemble[
        pred_columns
    ].mean(
        axis=1
    )


    ensemble[
        "predicted_age_sd"
    ] = ensemble[
        pred_columns
    ].std(
        axis=1,
        ddof=1
    )


    ensemble[
        "age_gap"
    ] = (

        ensemble[
            "predicted_age"
        ]

        -

        ensemble[
            "age"
        ]

    )


    # --------------------------------------------------------
    # Save ensemble
    # --------------------------------------------------------

    ensemble_path = os.path.join(

        output_dir,

        "ensemble_prediction.csv"

    )


    ensemble.to_csv(

        ensemble_path,

        index=False

    )


    # ========================================================
    # Overall metrics
    # ========================================================

    metrics = calculate_metrics(

        ensemble["age"].values,

        ensemble["predicted_age"].values

    )


    metrics["model"] = model_mode

    metrics["validation_type"] = (
        "right-trained_to_left-external"
    )


    metrics_df = pd.DataFrame(
        [metrics]
    )


    metrics_df.to_csv(

        os.path.join(

            output_dir,

            "external_metrics.csv"

        ),

        index=False

    )


    # --------------------------------------------------------
    # Save fold metrics based on prediction files
    # --------------------------------------------------------

    fold_metric_list = []


    for fold in FOLDS:

        col = (
            f"pred_fold{fold}"
        )


        fm = calculate_metrics(

            ensemble["age"].values,

            ensemble[col].values

        )


        fm["fold"] = fold

        fm["model"] = model_mode


        fold_metric_list.append(
            fm
        )


    pd.DataFrame(
        fold_metric_list
    ).to_csv(

        os.path.join(

            output_dir,

            "fold_metrics.csv"

        ),

        index=False

    )


    # ========================================================
    # Plots
    # ========================================================

    save_scatter_plot(

        ensemble,

        os.path.join(

            output_dir,

            "scatter_predicted_vs_chronological.png"

        ),

        (
            f"{model_mode}\n"
            "Right-trained → Left external validation"
        )

    )


    save_bland_altman(

        ensemble,

        os.path.join(

            output_dir,

            "bland_altman.png"

        ),

        (
            f"{model_mode}\n"
            "Right-trained → Left external validation"
        )

    )


    save_age_gap_plot(

        ensemble,

        os.path.join(

            output_dir,

            "age_gap_vs_age.png"

        ),

        (
            f"{model_mode}\n"
            "External validation"
        )

    )


    # ========================================================
    # Print final
    # ========================================================

    print()

    print(
        "=" * 90
    )

    print(
        f"{model_mode} "
        "EXTERNAL VALIDATION RESULT"
    )

    print(
        "=" * 90
    )


    print(
        f"N          : "
        f"{metrics['N']}"
    )

    print(
        f"MAE        : "
        f"{metrics['MAE']:.5f}"
    )

    print(
        f"RMSE       : "
        f"{metrics['RMSE']:.5f}"
    )

    print(
        f"R2         : "
        f"{metrics['R2']:.5f}"
    )

    print(
        f"Pearson r  : "
        f"{metrics['Pearson_r']:.5f}"
    )

    print()

    print(
        "Ensemble:"
    )

    print(
        ensemble_path
    )


    return metrics


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(

        description=(
            "External validation of "
            "right-trained ResNet101 "
            "models on left X-rays."
        )

    )


    parser.add_argument("--metadata_knee", type=str, default='data/external_left/single_left_knee_external.csv')
    parser.add_argument("--metadata_hip", type=str, default='data/external_left/single_left_hip_external.csv')
    parser.add_argument("--metadata_double", type=str, default='data/external_left/double_left_external.csv')


    parser.add_argument("--num_workers", type=int, default=DEFAULT_NUM_WORKERS)


    parser.add_argument(

        "--mode",

        type=str,

        default="all",

        choices=[
            "all",
            "knee",
            "hip",
            "double"
        ]

    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Device
    # --------------------------------------------------------

    if torch.cuda.is_available():

        device = torch.device(
            "cuda"
        )

        print(
            "Using:",
            torch.cuda.get_device_name(0)
        )

    else:

        device = torch.device(
            "cpu"
        )

        print(
            "Using CPU"
        )


    # --------------------------------------------------------
    # Create output
    # --------------------------------------------------------

    os.makedirs(

        EXTERNAL_DIR,

        exist_ok=True

    )


    # --------------------------------------------------------
    # Metadata mapping
    # --------------------------------------------------------

    tasks = {

        "knee": (

            "single_right_knee_51638",

            args.metadata_knee

        ),

        "hip": (

            "single_right_hip_51638",

            args.metadata_hip

        ),

        "double": (

            "double_right_51638",

            args.metadata_double

        )

    }


    # --------------------------------------------------------
    # Run
    # --------------------------------------------------------

    all_metrics = []


    for task_name, (
        model_mode,
        metadata_path
    ) in tasks.items():

        if (
            args.mode != "all"
            and
            args.mode != task_name
        ):

            continue


        metrics = run_external_validation(

            model_mode,

            metadata_path,

            device,

            args.num_workers

        )


        all_metrics.append(
            metrics
        )


    # ========================================================
    # Combined metrics
    # ========================================================

    if len(all_metrics) > 0:

        combined = pd.DataFrame(
            all_metrics
        )


        combined_path = os.path.join(

            EXTERNAL_DIR,

            "external_validation_summary.csv"

        )


        combined.to_csv(

            combined_path,

            index=False

        )


        print()

        print(
            "=" * 90
        )

        print(
            "ALL EXTERNAL VALIDATION RESULTS"
        )

        print(
            "=" * 90
        )

        print(
            combined[
                [
                    "model",
                    "N",
                    "MAE",
                    "RMSE",
                    "R2",
                    "Pearson_r"
                ]
            ].to_string(
                index=False
            )
        )

        print()

        print(
            "Summary:"
        )

        print(
            combined_path
        )


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    main()