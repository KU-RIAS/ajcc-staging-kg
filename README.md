# Large Language Model and Knowledge Graph–Driven AJCC Staging of Prostate Cancer Using Pathology Reports

## 📋 Overview
This repository implements an automated AJCC staging system that constructs a knowledge graph from extracted prostate cancer pathology report data and performs staging classification with consistency validation using Neo4j.


## 📁 Repository Structure
- `code/` - Source code
  - `kg_ajcc_staging.py`
- `data/` - Sample data
  - `sample_input.xlsx`

 
## 🛠️ Prerequisites

### 1. Python Dependencies
```bash
pip install pandas numpy neo4j openpyxl
```

### 2. Neo4j Database Setup

**Download and Install:**
- Download [Neo4j Desktop](https://neo4j.com/download/) 

**Create Database:**
1. Open Neo4j Desktop
2. Click **"Local instances"** → **"Create instance"**
3. Click **"Create"**

**Start Database:**
1. Click **"▶️Start"** on your created instance
2. **Before running the code:** Ensure your Neo4j credentials match those in `kg_ajcc_staging.py`:
```python
# Update these values if you used different settings
uri="bolt://localhost:7687"  # Connection URI
user="neo4j"                 # Username
password="password"          # Password
```

## 📊 Data Source
Sample data structure based on:

Grothey, B., & Tolkach, Y. (2024). Large Language Models for Extraction of Structured Data in Pathology (Version v3). Zenodo. https://doi.org/10.5281/zenodo.14175293
