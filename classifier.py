# =============================================================================
# classifier.py — Per-sample classification using SD and slope algorithms.
# =============================================================================

import config
from sd_algorithm import run_sd_algorithm
from slope_algorithm import run_slope_algorithm


def classify_sample(sample: dict) -> dict:
    """
    Run both detection algorithms on one sample and return combined results.

    Parameters
    ----------
    sample : dict with keys 'name', 'concentration', 'rgb_data'

    Returns
    -------
    dict with classification results for the sample
    """
    rgb_data = sample["rgb_data"]
    concentration = sample["concentration"]

    sd_result = run_sd_algorithm(rgb_data)
    slope_result = run_slope_algorithm(rgb_data)

    ground_truth = concentration >= config.UTI_POSITIVE_THRESHOLD_CFU

    final_rg = rgb_data["RG_ratio"][-1] if rgb_data.get("RG_ratio") else None
    final_rb = rgb_data["RB_ratio"][-1] if rgb_data.get("RB_ratio") else None
    final_gb = rgb_data["GB_ratio"][-1] if rgb_data.get("GB_ratio") else None

    return {
        "name": sample["name"],
        "concentration_cfu_ml": concentration,
        "ground_truth_positive": ground_truth,
        "sd_detected": sd_result["detected"],
        "sd_first_detection_min": sd_result["first_detection_min"],
        "sd_result": sd_result,
        "slope_detected": slope_result["detected"],
        "slope_first_detection_min": slope_result["first_detection_min"],
        "slope_result": slope_result,
        "final_RG_ratio": final_rg,
        "final_RB_ratio": final_rb,
        "final_GB_ratio": final_gb,
        "rgb_data": rgb_data,
    }


def classify_all(samples: list[dict]) -> list[dict]:
    """Classify all samples and return a list of result dicts."""
    results = []
    for sample in samples:
        result = classify_sample(sample)
        status = "POSITIVE" if result["sd_detected"] or result["slope_detected"] else "negative"
        print(f"  {sample['name']:30s}  SD: {result['sd_detected']}  "
              f"Slope: {result['slope_detected']}  [{status}]")
        results.append(result)
    return results
