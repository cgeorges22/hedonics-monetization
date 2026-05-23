# Hedonics v2.12

## Overview

Hedonics is a C++ simulation program that models firm dynamics in a macroeconomic setting. It simulates an economy with multiple firms and consumer markets, featuring product innovation, quality differentiation, and hedonic characteristics. v2.12 includes monetization mechanisms to investigate their welfare consequences. 

## Features

- **Hedonic Product Characterization**: Products are characterized by heterogeneous hedonic attributes 
- **Firm Innovation & R&D**: Endogenous firm product innovation dirven by endogenous research and development decisions
- **Monetization of Products**: by accumulation of artificial hedonic qualities
- **Consumer Optimization**: Search-based consumer behavior with utility-maximizing product selection
- **Market Segmentation**: Two classes of consumers and product markets with different dynamics

## Version History

- **v2.12** (May 2026): Added artificial and effective quality, monetization enhancements
- **v2.11** (June 2019): Macroeconomic gentrification experiments
- **v2.10** (June 2016): Single parameter changes during simulation runs
- **v2.8-v2.9**: Migration, market tracking, and consumer/firm segmentation
- **v2.0-v2.7**: Core features including labor accounting, complementarities, indivisibilities
- **v1.0-v1.9**: Initial hedonic product framework and R&D modeling

## Building

### Requirements
- C++ compiler (g++, clang, or MSVC)
- Standard C++ library

### Compilation

```bash
g++ -o hedonics hedonics-v2.12.cpp hedonics-utilities-v2.12.cpp -lm
```

Or with optimization flags:
```bash
g++ -O2 -o hedonics hedonics-v2.12.cpp hedonics-utilities-v2.12.cpp -lm
```

## Usage

```bash
./hedonics < input.txt
```

The program requires an `input.txt` file containing simulation parameters and configuration settings.

## Input Configuration

Parameters are specified in `input.txt` and control:
- Number of firms and consumers
- Market structures and segmentation
- R&D parameters and innovation rates
- Quality dynamics and hedonic attributes
- Output and labor specifications

## Output

The program generates detailed simulation results including:
- Firm profitability and market shares
- Quality evolution over time
- Consumer welfare metrics
- R&D investment patterns
- Market structure dynamics

## Authors

- **Chris Georges** (Primary developer, v2.0-v2.12)
- **Philip Ewing** (v1.9 parameter extraction)
- **Dane Johnson** (Original Firm Dynamics model, v1.0)

## Notes

This program is under active development. Features continue to be added and refined for economic modeling and monetization research.
