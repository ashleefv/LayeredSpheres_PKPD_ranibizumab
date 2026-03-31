import subprocess
import sys

VERSION_MAP = {
    "2024a": "24.1.4",
    "2024b": "24.2.2",
    "2025a": "25.1.2",
}

packages = [
    "torch==2.10.0",
    "gpytorch==1.15.1",
    "botorch==0.17.0",
    "pandas==2.3.3",
    "numpy==2.4.1",
    "scipy==1.16.3",
    "matplotlib==3.10.8",
]

for pkg in packages:
    print(f"Installing {pkg}...")
    subprocess.run([sys.executable, "-m", "pip", "install", pkg])

print("\nDetecting MATLAB version...")
result = subprocess.run(
    ["matlab", "-batch", "disp(version('-release'))"],
    capture_output=True, text=True
)
release = result.stdout.strip()

if release not in VERSION_MAP:
    print(f"MATLAB {release} not recognized. Install matlabengine manually.")
    sys.exit(1)

pkg_version = VERSION_MAP[release]
print(f"Detected MATLAB {release}, installing matlabengine=={pkg_version}...")
subprocess.run([sys.executable, "-m", "pip", "install", f"matlabengine=={pkg_version}"])
print("\nDone!")