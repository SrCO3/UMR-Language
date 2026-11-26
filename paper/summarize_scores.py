"""Recalculate the manuscript's paired pathway scores from its Appendix C table."""
import csv
from collections import Counter
from pathlib import Path


def summarize(rows):
    direct = sum(int(row["without_umr_full"]) for row in rows)
    umr = sum(int(row["with_umr_full"]) for row in rows)
    n = len(rows)
    return (
        f"{n:>3}   {direct:>3} ({100 * direct / n:>4.1f}%)"
        f"   {umr:>3} ({100 * umr / n:>4.1f}%)"
        f"   {umr - direct:>+3}   {100 * (umr - direct) / n:>+6.2f}"
    )


def main():
    path = Path(__file__).with_name("pathway-scores.csv")
    with path.open(encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))

    print("Set   N     Without UMR     With UMR    Gain   Delta (pp)")
    for name in ("A", "B", "C"):
        subset = [row for row in rows if row["case"].startswith(name)]
        print(f"{name:3} {summarize(subset)}")
    print(f"All {summarize(rows)}")

    pairs = Counter((row["without_umr_full"], row["with_umr_full"]) for row in rows)
    print(f"Both {pairs['1', '1']}, UMR only {pairs['0', '1']}, "
          f"without UMR only {pairs['1', '0']}, neither {pairs['0', '0']}.")
    required = sum(int(row["checkpoints"]) for row in rows)
    for label, field in (("Without UMR", "without_umr_matched"), ("With UMR", "with_umr_matched")):
        matched = sum(int(row[field]) for row in rows)
        print(f"{label} matched {matched}/{required} checkpoints.")
    statuses = Counter(row["with_umr_submission"] for row in rows)
    print(f"UMR submissions {statuses['submitted']}, no native submission "
          f"{statuses['no_native_submission']}, delivery unconfirmed {statuses['delivery_unconfirmed']}.")


if __name__ == "__main__":
    main()
