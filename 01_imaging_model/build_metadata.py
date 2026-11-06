# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 01_imaging_model
# 汇总影像清单、表型、协变量与结局，生成建模主表（替代早期 01/02/03 三步）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import re
import sys
import random
import hashlib

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold


# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT = "/path/to/project"

AGE_CSV = "/path/to/ref/Age_prediction/age_sex.csv"

IMAGE_ROOT = "/path/to/ref/Age_prediction"


IMAGE_DIRS = {

    "left_knee":
        os.path.join(
            IMAGE_ROOT,
            "Left_Ortho_Knee_normalized_fixedsize57647"
        ),

    "right_knee":
        os.path.join(
            IMAGE_ROOT,
            "Right_Ortho_Knee_mirror_normalized_fixedsize57497"
        ),

    "left_hip":
        os.path.join(
            IMAGE_ROOT,
            "Left_Femur_normalized_fixedsize61844"
        ),

    "right_hip":
        os.path.join(
            IMAGE_ROOT,
            "Right_Femur_mirror_normalized_fixedsize59927"
        )
}


OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "data",
    "model_metadata"
)


REPORT_DIR = os.path.join(
    PROJECT_ROOT,
    "results",
    "metadata_validation"
)


# ============================================================
# CV configuration
# ============================================================

N_SPLITS = 5

RANDOM_STATE = 42


# ============================================================
# Supported image extensions
# ============================================================

IMAGE_EXTENSIONS = {

    ".png",
    ".jpg",
    ".jpeg",
    ".bmp",
    ".tif",
    ".tiff"
}


# ============================================================
# Model definitions
# ============================================================

MODE_NAMES = [

    "single_left_knee",

    "single_right_knee",

    "single_left_hip",

    "single_right_hip",

    "double_knee",

    "double_hip",

    "double_left",

    "double_right",

    "four_joint"

]


# ============================================================
# Reproducibility
# ============================================================

random.seed(
    RANDOM_STATE
)

np.random.seed(
    RANDOM_STATE
)


# ============================================================
# Extract eid from filename
# ============================================================

def extract_eid(filename):

    """
    Example:

    4379817_20158_2_0_xxxxx_mirrored_final.png

    ->

    4379817
    """

    basename = os.path.basename(
        filename
    )

    match = re.match(
        r"^(\d+)_",
        basename
    )

    if match is None:

        return None

    return int(
        match.group(1)
    )


# ============================================================
# Build image index
# ============================================================

def build_image_index(
    image_dir,
    joint_name
):

    print()
    print("=" * 80)
    print(
        f"Scanning {joint_name}"
    )
    print(
        image_dir
    )
    print("=" * 80)


    if not os.path.isdir(
        image_dir
    ):

        raise FileNotFoundError(
            f"Image directory does not exist:\n"
            f"{image_dir}"
        )


    image_index = {}

    total_files = 0

    valid_files = 0

    invalid_files = 0


    for root, dirs, files in os.walk(
        image_dir
    ):

        for filename in files:

            total_files += 1


            ext = os.path.splitext(
                filename
            )[1].lower()


            if ext not in IMAGE_EXTENSIONS:

                continue


            eid = extract_eid(
                filename
            )


            if eid is None:

                invalid_files += 1

                continue


            full_path = os.path.join(
                root,
                filename
            )


            full_path = os.path.abspath(
                full_path
            )


            if not os.path.isfile(
                full_path
            ):

                continue


            valid_files += 1


            if eid not in image_index:

                image_index[eid] = []

            image_index[eid].append(
                full_path
            )


    # --------------------------------------------------------
    # Sort all duplicate images
    # --------------------------------------------------------

    for eid in image_index:

        image_index[eid] = sorted(
            image_index[eid]
        )


    duplicate_eids = {

        eid: paths

        for eid, paths
        in image_index.items()

        if len(paths) > 1

    }


    print(
        "Total files:",
        total_files
    )

    print(
        "Valid image files:",
        valid_files
    )

    print(
        "Invalid filename:",
        invalid_files
    )

    print(
        "Unique eid:",
        len(image_index)
    )

    print(
        "Eids with multiple images:",
        len(duplicate_eids)
    )


    return image_index


# ============================================================
# Select image
# ============================================================

def select_image(
    image_index,
    eid
):

    """
    If multiple images exist for one eid,
    use the first image after deterministic sorting.
    """

    paths = image_index.get(
        int(eid),
        []
    )


    if len(paths) == 0:

        return None


    return paths[0]


# ============================================================
# Check image
# ============================================================

def image_exists(path):

    if path is None:

        return False


    if not isinstance(
        path,
        str
    ):

        return False


    if path.strip() == "":

        return False


    return os.path.isfile(
        path
    )


# ============================================================
# Create participant fold
# ============================================================

def create_global_folds(
    df
):

    """
    Create ONE global participant-level 5-fold assignment.

    Important:
    The same eid keeps the same fold across all 9 models.

    This makes model comparison much cleaner.
    """


    work = df.copy()


    # --------------------------------------------------------
    # Age bins
    # --------------------------------------------------------

    n_bins = 10


    try:

        work["age_bin"] = pd.qcut(
            work["age"],
            q=n_bins,
            labels=False,
            duplicates="drop"
        )

    except Exception:

        work["age_bin"] = pd.cut(
            work["age"],
            bins=n_bins,
            labels=False
        )


    # --------------------------------------------------------
    # Sex + age bin
    # --------------------------------------------------------

    work["stratum"] = (

        work["sex"].astype(str)

        + "_"

        + work["age_bin"].astype(str)

    )


    # --------------------------------------------------------
    # Very small strata
    # --------------------------------------------------------

    stratum_counts = (
        work["stratum"]
        .value_counts()
    )


    small_strata = stratum_counts[
        stratum_counts < N_SPLITS
    ].index


    if len(small_strata) > 0:

        work.loc[
            work["stratum"].isin(
                small_strata
            ),
            "stratum"
        ] = "other"


    # --------------------------------------------------------
    # Stratified KFold
    # --------------------------------------------------------

    skf = StratifiedKFold(

        n_splits=N_SPLITS,

        shuffle=True,

        random_state=RANDOM_STATE

    )


    folds = np.zeros(
        len(work),
        dtype=int
    )


    X = np.zeros(
        len(work)
    )


    y = work[
        "stratum"
    ].values


    for fold_idx, (
        train_idx,
        val_idx
    ) in enumerate(
        skf.split(
            X,
            y
        ),
        start=1
    ):

        folds[val_idx] = fold_idx


    work["fold"] = folds.astype(
        int
    )


    return work[
        [
            "eid",
            "fold"
        ]
    ]


# ============================================================
# Build metadata row
# ============================================================

def create_row(
    row,
    fold,
    images
):

    return {

        "eid":
            int(row["eid"]),

        "age":
            float(row["age"]),

        "sex":
            int(row["sex"]),

        "fold":
            int(fold),

        **images

    }


# ============================================================
# Build all 9 metadata files
# ============================================================

def build_metadata(
    participants,
    image_indices
):

    results = {}

    reports = {}


    # ========================================================
    # Helper
    # ========================================================

    def get(
        joint,
        eid
    ):

        return select_image(
            image_indices[joint],
            eid
        )


    # ========================================================
    # SINGLE LEFT KNEE
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_knee",
            eid
        )

        right = get(
            "right_knee",
            eid
        )


        if image_exists(left):

            selected = left

            source = "left_knee"

        elif image_exists(right):

            selected = right

            source = "right_knee_fallback"

        else:

            reason[eid] = (
                "no_knee"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                selected

        })


    results[
        "single_left_knee"
    ] = pd.DataFrame(
        records
    )


    reports[
        "single_left_knee"
    ] = reason


    # ========================================================
    # SINGLE RIGHT KNEE
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_knee",
            eid
        )

        right = get(
            "right_knee",
            eid
        )


        if image_exists(right):

            selected = right

        elif image_exists(left):

            selected = left

        else:

            reason[eid] = (
                "no_knee"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                selected

        })


    results[
        "single_right_knee"
    ] = pd.DataFrame(
        records
    )


    reports[
        "single_right_knee"
    ] = reason


    # ========================================================
    # SINGLE LEFT HIP
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_hip",
            eid
        )

        right = get(
            "right_hip",
            eid
        )


        if image_exists(left):

            selected = left

        elif image_exists(right):

            selected = right

        else:

            reason[eid] = (
                "no_hip"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                selected

        })


    results[
        "single_left_hip"
    ] = pd.DataFrame(
        records
    )


    reports[
        "single_left_hip"
    ] = reason


    # ========================================================
    # SINGLE RIGHT HIP
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_hip",
            eid
        )

        right = get(
            "right_hip",
            eid
        )


        if image_exists(right):

            selected = right

        elif image_exists(left):

            selected = left

        else:

            reason[eid] = (
                "no_hip"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                selected

        })


    results[
        "single_right_hip"
    ] = pd.DataFrame(
        records
    )


    reports[
        "single_right_hip"
    ] = reason


    # ========================================================
    # DOUBLE KNEE
    #
    # left knee + right knee
    # both required
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_knee",
            eid
        )

        right = get(
            "right_knee",
            eid
        )


        if not image_exists(left):

            reason[eid] = (
                "missing_left_knee"
            )

            continue


        if not image_exists(right):

            reason[eid] = (
                "missing_right_knee"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                left,

            "image2":
                right

        })


    results[
        "double_knee"
    ] = pd.DataFrame(
        records
    )


    reports[
        "double_knee"
    ] = reason


    # ========================================================
    # DOUBLE HIP
    #
    # left hip + right hip
    # both required
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        left = get(
            "left_hip",
            eid
        )

        right = get(
            "right_hip",
            eid
        )


        if not image_exists(left):

            reason[eid] = (
                "missing_left_hip"
            )

            continue


        if not image_exists(right):

            reason[eid] = (
                "missing_right_hip"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                left,

            "image2":
                right

        })


    results[
        "double_hip"
    ] = pd.DataFrame(
        records
    )


    reports[
        "double_hip"
    ] = reason


    # ========================================================
    # DOUBLE LEFT
    #
    # Primary:
    # left knee + left hip
    #
    # Missing left knee -> right knee
    # Missing left hip  -> right hip
    #
    # Each joint is independently replaced.
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        lk = get(
            "left_knee",
            eid
        )

        rk = get(
            "right_knee",
            eid
        )

        lh = get(
            "left_hip",
            eid
        )

        rh = get(
            "right_hip",
            eid
        )


        # Knee
        if image_exists(lk):

            knee = lk

        elif image_exists(rk):

            knee = rk

        else:

            reason[eid] = (
                "no_knee"
            )

            continue


        # Hip
        if image_exists(lh):

            hip = lh

        elif image_exists(rh):

            hip = rh

        else:

            reason[eid] = (
                "no_hip"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                knee,

            "image2":
                hip

        })


    results[
        "double_left"
    ] = pd.DataFrame(
        records
    )


    reports[
        "double_left"
    ] = reason


    # ========================================================
    # DOUBLE RIGHT
    #
    # Primary:
    # right knee + right hip
    #
    # Missing right knee -> left knee
    # Missing right hip  -> left hip
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        lk = get(
            "left_knee",
            eid
        )

        rk = get(
            "right_knee",
            eid
        )

        lh = get(
            "left_hip",
            eid
        )

        rh = get(
            "right_hip",
            eid
        )


        # Knee
        if image_exists(rk):

            knee = rk

        elif image_exists(lk):

            knee = lk

        else:

            reason[eid] = (
                "no_knee"
            )

            continue


        # Hip
        if image_exists(rh):

            hip = rh

        elif image_exists(lh):

            hip = lh

        else:

            reason[eid] = (
                "no_hip"
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                knee,

            "image2":
                hip

        })


    results[
        "double_right"
    ] = pd.DataFrame(
        records
    )


    reports[
        "double_right"
    ] = reason


    # ========================================================
    # FOUR JOINT
    #
    # left knee
    # right knee
    # left hip
    # right hip
    #
    # ALL FOUR REQUIRED
    # ========================================================

    records = []

    reason = {}


    for _, row in participants.iterrows():

        eid = int(
            row["eid"]
        )

        lk = get(
            "left_knee",
            eid
        )

        rk = get(
            "right_knee",
            eid
        )

        lh = get(
            "left_hip",
            eid
        )

        rh = get(
            "right_hip",
            eid
        )


        missing = []


        if not image_exists(lk):

            missing.append(
                "left_knee"
            )


        if not image_exists(rk):

            missing.append(
                "right_knee"
            )


        if not image_exists(lh):

            missing.append(
                "left_hip"
            )


        if not image_exists(rh):

            missing.append(
                "right_hip"
            )


        if len(missing) > 0:

            reason[eid] = (
                "missing_"
                + "_".join(missing)
            )

            continue


        records.append({

            "eid":
                eid,

            "age":
                float(row["age"]),

            "sex":
                int(row["sex"]),

            "fold":
                int(row["fold"]),

            "image1":
                lk,

            "image2":
                rk,

            "image3":
                lh,

            "image4":
                rh

        })


    results[
        "four_joint"
    ] = pd.DataFrame(
        records
    )


    reports[
        "four_joint"
    ] = reason


    return results, reports


# ============================================================
# Validate metadata
# ============================================================

def validate_metadata(
    results
):

    print("\n")

    print("=" * 100)
    print("FINAL METADATA VALIDATION")
    print("=" * 100)


    summary = []


    for mode in MODE_NAMES:

        df = results[
            mode
        ]


        image_cols = [
            c
            for c in df.columns
            if c.startswith("image")
        ]


        invalid_rows = 0

        duplicate_eid = int(
            df["eid"]
            .duplicated()
            .sum()
        )


        for _, row in df.iterrows():

            for col in image_cols:

                if not image_exists(
                    row[col]
                ):

                    invalid_rows += 1

                    break


        fold_counts = (
            df["fold"]
            .value_counts()
            .sort_index()
        )


        print("\n")

        print(
            mode
        )

        print(
            "-" * 80
        )

        print(
            "Participants:",
            len(df)
        )

        print(
            "Unique eid:",
            df["eid"].nunique()
        )

        print(
            "Images per participant:",
            len(image_cols)
        )

        print(
            "Invalid image rows:",
            invalid_rows
        )

        print(
            "Duplicate eid:",
            duplicate_eid
        )

        print(
            "Fold counts:"
        )

        print(
            fold_counts.to_dict()
        )


        if invalid_rows != 0:

            raise RuntimeError(
                f"{mode} contains "
                f"invalid image paths."
            )


        if duplicate_eid != 0:

            raise RuntimeError(
                f"{mode} contains "
                f"duplicate eid."
            )


        summary.append({

            "mode":
                mode,

            "n":
                len(df),

            "unique_eid":
                df["eid"].nunique(),

            "n_images":
                len(image_cols),

            "invalid_image_rows":
                invalid_rows,

            "duplicate_eid":
                duplicate_eid,

            "fold1":
                int(
                    fold_counts.get(
                        1,
                        0
                    )
                ),

            "fold2":
                int(
                    fold_counts.get(
                        2,
                        0
                    )
                ),

            "fold3":
                int(
                    fold_counts.get(
                        3,
                        0
                    )
                ),

            "fold4":
                int(
                    fold_counts.get(
                        4,
                        0
                    )
                ),

            "fold5":
                int(
                    fold_counts.get(
                        5,
                        0
                    )
                )

        })


    return pd.DataFrame(
        summary
    )


# ============================================================
# Save exclusion reports
# ============================================================

def save_reports(
    participants,
    reports,
    results
):

    os.makedirs(
        REPORT_DIR,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary = []


    for mode in MODE_NAMES:

        result_df = results[
            mode
        ]

        excluded = reports[
            mode
        ]


        summary.append({

            "mode":
                mode,

            "age_available":
                len(participants),

            "included":
                len(result_df),

            "excluded":
                len(excluded),

            "excluded_percent":
                (
                    len(excluded)
                    / len(participants)
                    * 100
                )

        })


        if len(excluded) > 0:

            exclusion_df = pd.DataFrame({

                "eid":
                    list(
                        excluded.keys()
                    ),

                "reason":
                    list(
                        excluded.values()
                    )

            })


            exclusion_df.to_csv(

                os.path.join(

                    REPORT_DIR,

                    mode
                    + "_excluded.csv"

                ),

                index=False

            )


    summary_df = pd.DataFrame(
        summary
    )


    summary_df.to_csv(

        os.path.join(

            REPORT_DIR,

            "metadata_summary.csv"

        ),

        index=False

    )


    return summary_df


# ============================================================
# Main
# ============================================================

def main():
	
    print()

    print("=" * 100)

    print(
        "X-RAY AGE PREDICTION"
    )

    print(
        "Metadata generation"
    )

    print("=" * 100)


    # ========================================================
    # Check directories
    # ========================================================

    print("\nChecking image directories...")

    for joint, path in IMAGE_DIRS.items():

        print(
            f"{joint:15s}: "
            f"{path}"
        )


        if not os.path.isdir(path):

            raise FileNotFoundError(
                f"Directory does not exist:\n"
                f"{path}"
            )


    # ========================================================
    # Load age CSV
    # ========================================================

    print("\n")

    print(
        "Loading:",
        AGE_CSV
    )


    age_df = pd.read_csv(
        AGE_CSV
    )


    required_cols = [
        "eid",
        "sex",
        "age_i2"
    ]


    for col in required_cols:

        if col not in age_df.columns:

            raise ValueError(
                f"Missing column: {col}"
            )


    print(
        "Original rows:",
        len(age_df)
    )


    # ========================================================
    # Clean eid
    # ========================================================

    age_df["eid"] = pd.to_numeric(
        age_df["eid"],
        errors="coerce"
    )


    age_df["age_i2"] = pd.to_numeric(
        age_df["age_i2"],
        errors="coerce"
    )


    age_df["sex"] = pd.to_numeric(
        age_df["sex"],
        errors="coerce"
    )


    # Remove invalid eid
    age_df = age_df[
        age_df["eid"].notna()
    ].copy()


    age_df["eid"] = (
        age_df["eid"]
        .astype(np.int64)
    )


    # ========================================================
    # Keep participants with age
    # ========================================================

    participants = age_df[
        age_df["age_i2"].notna()
    ].copy()


    participants = participants[
        participants["sex"].notna()
    ].copy()


    participants = participants.rename(

        columns={

            "age_i2":
                "age"

        }

    )


    # --------------------------------------------------------
    # Duplicate eid check
    # --------------------------------------------------------

    if participants[
        "eid"
    ].duplicated().any():

        duplicated = participants[
            participants[
                "eid"
            ].duplicated(
                keep=False
            )
        ]


        print(
            duplicated
        )


        raise RuntimeError(
            "Duplicate eid in age_sex.csv"
        )


    print(
        "Participants with valid age:",
        len(participants)
    )


    # ========================================================
    # Global 5-fold
    # ========================================================

    print("\n")

    print(
        "Creating global "
        "participant-level 5-fold..."
    )


    fold_df = create_global_folds(
        participants
    )


    participants = participants.merge(

        fold_df,

        on="eid",

        how="left",

        validate="one_to_one"

    )


    if participants[
        "fold"
    ].isna().any():

        raise RuntimeError(
            "Some participants "
            "do not have a fold."
        )


    participants["fold"] = (
        participants["fold"]
        .astype(int)
    )


    print(
        "Global fold counts:"
    )

    print(
        participants[
            "fold"
        ].value_counts()
        .sort_index()
        .to_dict()
    )


    # ========================================================
    # Build image indices
    # ========================================================

    image_indices = {}


    for joint, path in IMAGE_DIRS.items():

        image_indices[joint] = (
            build_image_index(
                path,
                joint
            )
        )


    # ========================================================
    # Image availability summary
    # ========================================================

    print("\n")

    print("=" * 100)

    print(
        "IMAGE AVAILABILITY"
    )

    print("=" * 100)


    availability = []


    for joint in [
        "left_knee",
        "right_knee",
        "left_hip",
        "right_hip"
    ]:

        index = image_indices[
            joint
        ]


        n = sum(

            int(
                eid in index
            )

            for eid
            in participants["eid"]

        )


        availability.append({

            "joint":
                joint,

            "participants_with_image":
                n,

            "participants_without_image":
                len(participants) - n

        })


        print(
            f"{joint:15s}: "
            f"{n:8d} "
            f"available / "
            f"{len(participants)-n:8d} "
            f"missing"
        )


    pd.DataFrame(
        availability
    ).to_csv(

        os.path.join(

            REPORT_DIR,

            "joint_availability.csv"

        ),

        index=False

    )


    # ========================================================
    # Build metadata
    # ========================================================

    print("\n")

    print(
        "Building 9 metadata files..."
    )


    results, reports = build_metadata(

        participants,

        image_indices

    )


    # ========================================================
    # Save metadata
    # ========================================================

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True
    )


    for mode in MODE_NAMES:

        df = results[
            mode
        ]


        output_path = os.path.join(

            OUTPUT_DIR,

            mode + ".csv"

        )


        df.to_csv(

            output_path,

            index=False

        )


        print()

        print(
            f"{mode:25s}"
            f"N = {len(df):8d}"
        )

        print(
            output_path
        )


    # ========================================================
    # Save exclusion reports
    # ========================================================

    summary_df = save_reports(

        participants,

        reports,

        results

    )


    # ========================================================
    # Validate
    # ========================================================

    validation_df = validate_metadata(
        results
    )


    validation_df.to_csv(

        os.path.join(

            REPORT_DIR,

            "metadata_validation.csv"

        ),

        index=False

    )


    # ========================================================
    # Final summary
    # ========================================================

    print("\n")

    print("=" * 100)

    print(
        "FINAL SUMMARY"
    )

    print("=" * 100)


    print(
        summary_df.to_string(
            index=False
        )
    )


    print("\n")

    print(
        "Metadata directory:"
    )

    print(
        OUTPUT_DIR
    )


    print("\n")

    print(
        "Validation reports:"
    )

    print(
        REPORT_DIR
    )


    print("\n")

    print("=" * 100)

    print(
        "METADATA GENERATION COMPLETE"
    )

    print("=" * 100)
    # ========================================================
    # Create output directories
    # ========================================================
	
    os.makedirs(OUTPUT_DIR,exist_ok=True)
	
    os.makedirs(REPORT_DIR,exist_ok=True)

# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    main()