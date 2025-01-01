# Analysis of Emergency Room (ER) Visits in California (2011–2019) 🏥

## Description
This project presents an in-depth analysis of Emergency Room (ER) visit trends across California from 2011 to 2019. Using **Bayesian hierarchical models** implemented with **R-INLA**, the study examines fluctuations in ER visitation rates while accounting for spatial and spatio-temporal dynamics.

A key focus is on the impact of environmental and demographic factors, such as:
- **Atmospheric PM2.5 concentration**: A known health hazard correlated with respiratory and cardiovascular health issues.
- **Demographic composition**: Specifically, the proportion of the white population in each county and its influence on healthcare utilization patterns.

---

## Features
- **Bayesian Hierarchical Modeling**: Utilized R-INLA for robust statistical modeling.
- **Spatio-Temporal Analysis**: Explored regional variability and temporal trends in ER visitation rates.
- **Impact Assessment**:
  - **Environmental Factors**: Correlated PM2.5 levels with ER visit counts.
  - **Demographics**: Investigated the relationship between county demographics and healthcare utilization.
- **Policy-Oriented Insights**: Identified trends and determinants to inform public health interventions and emergency medical service optimization.

---

## My Contributions
1. **Model Development**: Designed and implemented Bayesian hierarchical models using R-INLA to incorporate spatial and temporal dynamics.
2. **Variable Integration**: Integrated PM2.5 concentrations and demographic proportions into the models for a comprehensive analysis.
3. **Result Interpretation**: Identified significant temporal trends and evaluated the added value of spatial components in model performance.
4. **Findings Communication**: Delivered actionable insights for policymakers, healthcare providers, and public health professionals.

---

## Tech Stack
- **Programming Language**: R
- **Key Libraries**:
  - `INLA` for Bayesian hierarchical modeling
  - `sp` and `sf` for spatial data manipulation
  - `ggplot2` for visualization
  - `dplyr` for data processing
- **Tools**: RStudio, Markdown

---

## Key Findings
- **Environmental Impact**: PM2.5 levels showed a significant correlation with ER visit counts, emphasizing the need for environmental health interventions.
- **Demographic Influence**: The demographic composition, particularly the proportion of the white population, influenced healthcare utilization patterns.
- **Spatio-Temporal Modeling**:
  - Spatial components alone did not significantly improve model performance.
  - Temporal dynamics highlighted counties with significant changes in ER visit trends over time.

---

