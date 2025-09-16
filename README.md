# Meta-Analysis: Maternal Toxicant Effects on Fish Reproductive Success

## Overview

This repository contains the complete dataset and analysis files for a meta-analysis examining how mercury (Hg) and polychlorinated biphenyls (PCBs) in maternal fish affect toxicant transfer to eggs and subsequent offspring success.

## Research Questions

1. **Maternal Size Effects**: How does maternal fish size relate to egg toxicant concentrations?
2. **Maternal Muscle Transfer**: How do maternal muscle toxicant levels relate to egg concentrations?  
3. **Offspring Success**: How do egg toxicant levels affect hatching success, mortality, and growth?

## Dataset Description

- **Studies Included**: 42 peer-reviewed studies (1993-2020)
- **Toxicants**: Mercury (Hg) and Polychlorinated Biphenyls (PCBs)
- **Endpoints**: Size, Muscle, Hatching, Mortality, Growth

## Methodology

### Search Strategy
Systematic literature searches conducted via ISI Web of Science:
- **Maternal Size Analysis** (11/23/2020): 713 papers + 413 forward citations reviewed, 9 used
- **Maternal Muscle Analysis** (12/29/2020): 571 papers reviewed, 19 used  
- **Offspring Success Analysis** (3/9/2021): 441 papers reviewed, 21 used

### Two Analytical Approaches

**Conservative Approach**:
- One effect size per study per analysis
- Largest sample size selected for multiple populations
- Random selection for equal sample sizes
- Single treatment level for experimental studies

**Inclusive Approach**:
- Multiple effect sizes per study allowed
- Separate effect size for each population/site
- One effect size per population per contaminant
- Multiple treatment levels included

## File Structure

```
├── Data/
│   ├── VBAProject.xlsx           # Complete dataset with all analyses
│   ├── Studies_Table/            # Master bibliography and study characteristics
│   ├── Search_Documentation/     # Literature search terms and results
│   ├── Maternal_Size_Analysis/   # Size vs egg toxicant relationships
│   ├── Maternal_Muscle_Analysis/ # Muscle vs egg toxicant relationships
│   ├── Hatching_Analysis/        # Egg toxicants vs hatching success
│   ├── Mortality_Analysis/       # Egg toxicants vs offspring mortality
│   └── Growth_Analysis/          # Egg toxicants vs offspring growth
├── Results/
│   ├── Statistical_Summaries/    # Meta-analysis results and effect sizes
│   └── Publication_Bias_Tests/   # Egger's tests and funnel plots
└── Documentation/
    ├── README.md                 # This file
    ├── Metadata.md              # Detailed dataset documentation
    └── Methods.md               # Statistical methods and calculations
```

## Data Types

### Study Characteristics
- **Maternal Size Studies**: Field collections measuring maternal length/age vs egg toxicants
- **Maternal Muscle Studies**: Field collections (+ 2 lab studies) measuring muscle vs egg toxicants  
- **Offspring Success Studies**: Laboratory experiments with controlled toxicant exposures

### Effect Size Calculations
- **Correlational analyses**: Fisher's z-transformed correlations for maternal effects
- **Experimental analyses**: Hedge's d for mean differences between treatment/control groups
- **Conversions**: Methods included for converting between correlation and mean difference effect sizes

## Key Variables

- **Toxicant Concentrations**: Wet weight basis (µg/kg ww), dry weights converted using 80% water content
- **Mercury**: Total mercury and methylmercury treated as combined category
- **PCBs**: All PCB congeners grouped, total PCBs preferred over individual congeners
- **Fish Size**: Length preferred over mass or age when available
- **Sample Units**: Individual fish for maternal analyses, group replicates for offspring studies

## Statistical Approach

- **Weighted mean effect sizes**: Inverse variance weighting
- **Confidence intervals**: 95% CI using Student's t-distribution
- **Publication bias**: Egger's tests and funnel plot analysis
- **Heterogeneity**: Variance estimates and bias correction methods

## Usage Notes

### Data Access
The main dataset is contained in `VBAProject.xlsx` with separate sheets for each analysis type and approach (conservative vs inclusive).

### Replication
All effect size calculations follow Borenstein et al. (2009) and Gurevitch et al. (2001) methods. Statistical formulas and conversion factors are documented in the Methods section.

### Limitations
- Freshwater fish bias (limited marine species)
- Laboratory exposure levels may exceed environmental concentrations
- Potential publication bias in offspring success analyses

## Citation

If using this dataset, please cite:

Sosa, B.M., Malinowski, C.R., & Höök, T.O. (2024). Implications of Maternal Toxicant Effects on Size-Dependent Fisheries' Management: A Meta-Analysis. [Journal Information]

## Authors

- **Brandon M. Sosa** - Florida International University, Purdue University
- **Christopher R. Malinowski** - Purdue University, Ocean First Institute  
- **Tomas O. Höök** - Purdue University, Illinois-Indiana Sea Grant


## Contact

For questions about the dataset or methodology, please contact:
- Brandon M. Sosa: [bsosa012@fiu.edu]


## Acknowledgments

This research was supported by Purdue University's Department of Forestry and Natural Resources, the Illinois-Indiana Sea Grant College Program, and the US Department of Agriculture (Hatch program).
