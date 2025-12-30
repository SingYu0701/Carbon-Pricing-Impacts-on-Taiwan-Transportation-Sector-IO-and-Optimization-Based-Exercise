# Carbon-Pricing-Impacts-on-Taiwan-Transportation-Sector-IO-and-Optimization-Based-Exercise
Final report of Interindustry Economics, Dec 2025 @ NCKU Resource Engineering 

![Made with R](https://img.shields.io/badge/Made%20with-R-276DC3?logo=r&logoColor=white)

## Abstract

This repository presents an academic exercise applying input–output (I–O) economics, carbon cost accounting, and linear programming optimization to analyze the potential economic impacts of carbon pricing on Taiwan’s transportation-related industries. The study focuses on how carbon costs affect sectoral output, propagate through inter-industry linkages, and are absorbed by different components of value added under alternative policy objectives.

The project was conducted as part of a graduate-level course on industrial linkage and input–output analysis. Its primary purpose is methodological exploration and analytical training, rather than empirical policy evaluation. All results should therefore be interpreted as model-based simulations under stylized assumptions.

## 1. Research Motivation

Transportation-related industries play a critical role in modern economies due to their strong upstream and downstream linkages, high energy intensity, and strategic importance for logistics, tourism, and international trade. Among these activities, aviation and freight transportation are particularly exposed to carbon pricing policies because of their relatively high emissions per unit of output.

In Taiwan, however, official input–output tables aggregate aviation, land transport, maritime transport, and logistics services into a single transportation sector, limiting the feasibility of disaggregated carbon impact analysis. This project therefore adopts an aggregated transportation sector as the core unit of analysis, aiming to examine:

- **The magnitude of direct carbon costs under alternative carbon price scenarios**

- **The output and demand responses induced by cost pass-through and demand elasticity**

- **The inter-industry transmission mechanisms of carbon cost shocks**

- **The distribution of absorbed carbon costs across labor income, capital consumption, and operating surplus**

Despite data aggregation constraints, the framework provides a structured and transparent approach for understanding the systemic role of transportation industries in a carbon-priced economy.

## 2. Analytical Framework
### 2.1 Input–Output Model

Input–Output Analysis is a core tool for measuring the interdependence among industries (Miller & Blair, 2009). In this study, the original 35 sectors are aggregated into 12 sectors. Accordingly, we use the intermediate input matrix aggregated to 12 sectors, which can be expressed as the sum of the original 35×35 matrix:

$$
Z_{12x12}[m,n] = \sum_{i \in S_m} \sum_{j \in S_n} Z_{35}[i,j], \quad m,n = 1,\dots,12
$$

where $$S_m$$ and $$S_n$$ denote the sets of original 35 subsectors corresponding to the m-th and n-th aggregated sectors, respectively. Final demand and carbon emission data are summed in the same manner to obtain 12-sector vectors. Let $$x∈R^{12}$$ denote the total output vector of each sector. The technical coefficient matrix $$A$$ is then defined as:

$$
A = Z X^{-1}
$$

Leontief (1986) proposed the Leontief inverse matrix $$L$$ to calculate the sensitivity of output to final demand, capturing both the direct effect of increased sector demand on its own output and the indirect effect via multi-round upstream purchases of intermediate inputs. This allows quantifying the economic ripple effects caused by changes in final demand due to rising carbon prices. The Leontief inverse is defined as:

$$L=(I−A)^{−1}$$

The change in output $$Δx$$ in response to a shock in final demand $$Δf$$ is then calculated as:

$$Δx=LΔf$$

This framework enables the simulation of **multiplier effects**, capturing how shocks to one sector propagate through the production network.

### 2.2 Carbon Cost Estimation and Demand Adjustment

Sectoral carbon costs are calculated by combining sector-level CO₂ emissions with assumed carbon prices. Two behavioral responses are modeled:

$$
C_{i}^{direct} = p_{CO_2} \bullet E_i
$$

​where $$E_i$$ denotes CO₂ emissions of sector $$i$$, and $$p_{CO_2}$$ is the carbon price (NTD per million tons CO₂).

- Cost absorption by firms, which compresses profits or factor incomes

$$
C_{i}^{absorb} =(1-τ)C_{i}^{direct}
$$

- Cost pass-through to prices, which affects final demand depending on assumed price elasticity

$$
C_{i}^{passed} = τC_{i}^{direct}
$$

The proportions of cost pass-through and profit absorption are $$τ$$ and $$1−τ$$.

Considering the price elasticity of final demand $$ε$$, the change in final demand can be approximated as:

$$Δf_i^{demand}=f_i[(1+price_increasei)^ε−1], price_i^{increase}=C_i^{passed}/x_i$$

This step captures the reduction in consumer demand due to carbon pricing and its **ripple effects on the supply chain via the Leontief model**.

### 2.3 Linear Programming Optimization

To explore policy-relevant trade-offs, a linear programming (LP) model is implemented under multiple objective functions, including:

- **Minimization of absorbed carbon costs**

- **Minimization of total economic costs**

- **Maximization of total output**

- **Minimization of labor income losses**

Output constraints and minimum demand requirements are imposed to ensure economically feasible solutions. This component highlights how different policy priorities may lead to distinct sectoral adjustment paths.

## 3. Empirical Results and Model-Based Findings
### 3.1 Direct Carbon Cost Distribution

Simulation results indicate that **energy-intensive and transportation-related** sectors bear the highest direct carbon costs under both low and high carbon price scenarios. In particular:

- Transportation emerges as one of the top contributors to total carbon costs, reflecting its high energy consumption and extensive upstream dependencies.

- Energy supply and basic materials sectors exhibit even higher absolute carbon costs, emphasizing their upstream exposure within the production network.

<img width="2560" height="1317" alt="圖片" src="https://github.com/user-attachments/assets/a88ff834-79be-4f60-8a0a-e1e606457de0" />



These findings confirm the structural vulnerability of transportation and heavy industry to carbon pricing mechanisms.

### 3.2 Output Responses under Cost Pass-Through Scenarios

When carbon costs are fully or partially passed through to prices, the model produces positive output-equivalent adjustments in several sectors. This result does not imply increased real consumption but instead reflects the equivalent revenue required to maintain production levels under higher unit costs.

Key patterns include:

- Higher carbon prices combined with full cost pass-through generate the largest output-equivalent effects **due to multiplier amplification**.

- **The transportation sector exhibits strong indirect effects on upstream industries, particularly energy supply and materials manufacturing**.

- Output responses are highly sensitive to assumptions regarding demand elasticity, **underscoring the importance of consumer behavior** in carbon pricing outcomes.

<img width="1072" height="766" alt="圖片" src="https://github.com/user-attachments/assets/76e9619d-218b-497a-886d-4c42a3b3ee66" />



### 3.3 Effects of Demand Elasticity

Incorporating price elasticity of demand moderates the magnitude of output-equivalent adjustments. While most sectors experience reduced output under elastic demand assumptions, certain **upstream sectors (e.g., energy and basic materials)** continue to display positive output responses.

<img width="2560" height="1317" alt="圖片" src="https://github.com/user-attachments/assets/1ac4db05-57af-4544-89fc-c707a0612d15" />


This counterintuitive result arises from **strong inter-industry multiplier effects**, whereby upstream input demand remains robust despite contraction in final consumption.

### 3.4 Value-Added Absorption of Carbon Costs

Absorbed carbon costs are allocated across three value-added components: labor compensation, capital consumption, and operating surplus. Results reveal pronounced heterogeneity:

- Profit compression dominates in energy- and capital-intensive sectors.

- Labor income absorbs a larger share of carbon costs in labor-intensive industries, raising concerns about wage pressure.

- Capital-intensive sectors experience higher burdens through increased capital consumption, potentially affecting long-term investment incentives.

<img width="2560" height="1317" alt="圖片" src="https://github.com/user-attachments/assets/9cc1cd06-6f77-47ee-8fac-3826ffabfd9e" />



These patterns illustrate how carbon pricing can influence income distribution within and across sectors.

### 3.5 Inter-Industry Transmission of Transportation Sector Shocks

A targeted simulation of a unit demand shock in the transportation sector demonstrates its role as a systemic transmission hub. Output responses are particularly pronounced in:

- Energy supply sectors

- Basic materials manufacturing

- General manufacturing and service industries

<img width="2560" height="1317" alt="圖片" src="https://github.com/user-attachments/assets/fd58fb42-a462-4ac2-9882-41e4fcee0bae" />


This result highlights the centrality of transportation in sustaining production networks and the potential for **amplified system-wide impacts arising from sector-specific cost shocks**.

### 3.6 Optimization Outcomes under Alternative Policy Objectives

Linear programming simulations reveal that different policy objectives lead to markedly different sectoral adjustment patterns:

1.Objectives focused on minimizing absorbed costs favor output reductions in energy-intensive sectors.

2.Output-maximization objectives maintain production levels but shift cost burdens toward labor or capital.

3.Labor-protection objectives reallocate carbon costs away from labor income, often at the expense of higher capital or profit losses.

<img width="2560" height="1317" alt="圖片" src="https://github.com/user-attachments/assets/4eac9640-b886-4bff-8f59-2c485f6f2077" />


These results underscore the importance of explicitly defining policy priorities when designing carbon pricing and compensation mechanisms.

## 4. Discussion

The findings demonstrate that the economic impacts of carbon pricing are **highly nonlinear and structurally mediated by inter-industry linkages and factor composition**. Transportation plays a pivotal role in amplifying cost shocks, while upstream sectors exhibit strong multiplier-driven resilience.

From a methodological perspective, the integration of input–output analysis with optimization techniques provides a flexible framework for exploring trade-offs among competing policy objectives.

## 5. Limitations

- Sector aggregation conceals heterogeneity within transportation and manufacturing subsectors

- Behavioral parameters are assumed rather than empirically estimated

- The model is static and does not capture dynamic adjustment or technological change

These limitations are consistent with the educational scope of the project.

## 6. Intended Use and Disclaimer

This repository is intended for **academic demonstration and methodological learning**. All results are scenario-based simulations and should not be interpreted as empirical forecasts or policy recommendations.

## 7. Reference

  Dantzig, G. B. (1951). Maximization of a linear function of variables subject to linear inequalities. Activity Analysis of Production and Allocation, 13, 339–347.
  
  EU Emissions Trading System (EU ETS) reports
  
  Leontief, W. (1986). Input-Output Economics. Oxford University Press.
  
  Ma, N., Li, H., Zhang, J., Han, X., Feng, S., & Arif, A. (2021). The short-term price effects and transmission mechanism of CO₂ cost pass‑through in China: A partial transmission model. Resources Policy, 70, 101972. https://doi.org/10.1016/j.resourpol.2020.101972
  
  Miller, R. E., & Blair, P. D. (2009). Input-Output Analysis: Foundations and Extensions. Cambridge University Press.
  
  Zhang, H., Hewings, G. J. D., & Zheng, X. (2019). The effects of carbon taxation in China: An analysis based on energy input-output model in hybrid units. Energy Policy, 128, 223–234. https://doi.org/10.1016/j.enpol.2018.12.045
