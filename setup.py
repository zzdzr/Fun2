from setuptools import setup, find_packages
from pathlib import Path

setup(
    name="fun2",
    version="0.1.0",
    description="Reinforcement-learning framework to track replication fountain structures.",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=Path("./env/requirements.txt").read_text().splitlines(),
    entry_points={
        "console_scripts": ["fun2=fun2.cli.implement_fun2:main"],
    },
)
