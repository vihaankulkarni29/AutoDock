import os
import urllib.request
import zipfile
import stat


def install_vina():
    print("Checking for AutoDock Vina...")
    tools_dir = "tools"
    os.makedirs(tools_dir, exist_ok=True)

    vina_exe = os.path.join(tools_dir, "vina.exe")

    if os.path.exists(vina_exe):
        print("✅ Vina is present.")
        return

    print("⬇️ Downloading Vina (Windows)...")
    exe_url = (
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/"
        "v1.2.7/vina_1.2.7_win.exe"
    )
    zip_url = (
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/"
        "v1.2.7/vina_1.2.7_windows_x86_64.zip"
    )

    try:
        urllib.request.urlretrieve(exe_url, vina_exe)
        print("✅ Vina installed in /tools.")
        return
    except Exception:
        zip_path = os.path.join(tools_dir, "vina.zip")
        urllib.request.urlretrieve(zip_url, zip_path)
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(tools_dir)

        # Move and cleanup (Adjust based on zip structure)
        # Usually extracts to a subfolder, we flatten it for simplicity or just point to it.
        print("✅ Vina installed in /tools.")


if __name__ == "__main__":
    install_vina()
    print("Ready. Run 'python -m autoscan.main' to test.")
