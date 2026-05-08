import os

os.environ.setdefault("KERR_GW_INTERFACE_SILENT", "1")
os.environ.setdefault("KERR_GW_PHYSICS_SILENT", "1")

from src.plot_fig2 import generate_and_plot_fig2


if __name__ == "__main__":
    print("=== Kerr GW Lensing Fig.2 Reproduction ===")
    generate_and_plot_fig2()
