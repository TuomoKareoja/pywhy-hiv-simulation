# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Description
This repository is for me to learn how to use the PyWhy causal framework. I want to test different methods and their limits by applying them to simulated
data created with whynot library that can create 'realistic' synthetic data for specific domains. The analyses with PyWhy are all in notebooks to make it easy to read the analyses directly from GitHub

## Environment
- Use conda environment from environment.yml: `conda activate pywhy-hiv-simulation`

## Code Style
- Follow PEP 8 conventions
- Use type annotations
- Organize imports: stdlib, third-party, local
- Naming: snake_case for functions/variables, CamelCase for classes
- Docstrings: NumPy style documentation
- Error handling: Use explicit exception handling with appropriate error messages
- Maximum line length: 120 characters (Black formatter limit increased with config)
