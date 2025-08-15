#!/usr/bin/env python3
"""
install_requirements.py
------------------------

A helper script to install all Python dependencies listed in:
Requirements/requirements.txt

Usage:
    python install_requirements.py
"""

import subprocess
import sys
import os

def main():
    requirements_path = os.path.join("Requirements", "requirements.txt")
    
    if not os.path.exists(requirements_path):
        sys.exit(f"❌ Requirements file not found at: {requirements_path}")
    
    print(f"📦 Installing packages from {requirements_path}...\n")
    
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", requirements_path])
        print("\n✅ All requirements installed successfully!")
    except subprocess.CalledProcessError:
        sys.exit("\n❌ Installation failed. Please check your internet connection or the package versions.")

if __name__ == "__main__":
    main()
