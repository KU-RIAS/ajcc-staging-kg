# Knowledge Graph-based AJCC Staging System for Prostate Cancer Pathology Reports

import pandas as pd
import numpy as np
from neo4j import GraphDatabase
import warnings
import sys
import datetime
import os
warnings.filterwarnings('ignore')



class TeeOutput:
    """
    Console output redirector that writes to both console and file
    """
    def __init__(self, log_file_path):
        self.terminal = sys.stdout
        self.log_file = open(log_file_path, 'w', encoding='utf-8')
        self.log_file.write(f"=== Prostate Cancer Knowledge Graph Log ===\n")
        self.log_file.write(f"Start Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.log_file.write("=" * 50 + "\n\n")
        
    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()
        
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()
        
    def close(self):
        self.log_file.write(f"\n" + "=" * 50 + "\n")
        self.log_file.write(f"End Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.log_file.write("=== Log End ===\n")
        self.log_file.close()


def setup_logging(log_file_path):
    """Setup logging to both console and file"""
    tee = TeeOutput(log_file_path)
    sys.stdout = tee
    return tee


class ProstateCancerKnowledgeGraph:
    """
    Prostate Cancer Knowledge Graph Builder with Single-Pass AJCC Classification
    
    Features:
    - Single-pass classification using available data only
    - Comprehensive consistency validation
    - Neo4j knowledge graph construction
    - Detailed classification failure tracking
    - Consistency error reporting and Excel export
    """

    def __init__(self, uri="bolt://localhost:7687", user="neo4j", password="password"):
        """Initialize Neo4j connection"""
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        """Close Neo4j connection"""
        self.driver.close()

    # Database Setup
    
    def clear_database(self):
        """Clear all nodes and relationships from database"""
        with self.driver.session() as session:
            session.run("MATCH (n) DETACH DELETE n")

    def create_constraints(self):
        """Create unique constraints for nodes"""
        with self.driver.session() as session:
            constraints = [
                "CREATE CONSTRAINT patient_id IF NOT EXISTS FOR (p:Patient) REQUIRE p.id IS UNIQUE",
                "CREATE CONSTRAINT report_id  IF NOT EXISTS FOR (r:PathologyReport) REQUIRE r.id IS UNIQUE",
                "CREATE CONSTRAINT stage_code IF NOT EXISTS FOR (s:TNM_Stage) REQUIRE s.code IS UNIQUE",
                "CREATE CONSTRAINT ajcc_stage IF NOT EXISTS FOR (a:AJCC_Stage) REQUIRE a.stage IS UNIQUE",
            ]
            for constraint in constraints:
                try:
                    session.run(constraint)
                except Exception:
                    pass  # Constraint already exists

    # Data Processing Utilities
    
    @staticmethod
    def _psa_to_float(val):
        """
        Convert PSA value to float
        
        Args:
            val: PSA value (string, int, float, or None)
            
        Returns:
            float or None: Converted PSA value, None if conversion fails
        """
        if pd.isna(val):
            return None
        if isinstance(val, (int, float)):
            return float(val)
        try:
            return float(str(val).strip().replace(",", ""))
        except Exception:
            return None

    @staticmethod
    def _is_missing(value) -> bool:
        """
        Check if a value is missing or not mentioned
        """
        if pd.isna(value):
            return True
        if isinstance(value, str):
            cleaned = value.strip().lower()
            if cleaned in {"not mentioned", "-", "", "nan", "null", "unknown", "none"}:
                return True
        return False

    def calculate_ajcc_stage(self, row):
        """
        Single-pass AJCC Stage calculation.
        If unknown, include the reason in parentheses.
        """
        error_reasons = []
        
        try:
            # Grade Group parsing
            gg_raw = str(row["WHO Grade Group"]).strip().lower()
            gg = None
            
            if pd.isna(row["WHO Grade Group"]) or gg_raw in ['nan', '']:
                error_reasons.append("Missing Grade Group")
            elif "gg" in gg_raw:
                try:
                    gg = int(gg_raw.replace("gg", "").strip())
                    if not (1 <= gg <= 5):
                        error_reasons.append("Invalid Grade Group range")
                        gg = None
                except:
                    error_reasons.append("Invalid Grade Group format")
            elif gg_raw.isdigit() and 1 <= int(gg_raw) <= 5:
                gg = int(gg_raw)
            else:
                error_reasons.append("Invalid Grade Group")
            
            # PSA parsing
            psa_val = None
            if pd.isna(row["Serum PSA"]):
                error_reasons.append("Missing PSA")
            else:
                try:
                    psa_val = self._psa_to_float(row["Serum PSA"])
                    if psa_val is None or psa_val < 0:
                        error_reasons.append("Invalid PSA value")
                        psa_val = None
                except:
                    error_reasons.append("Invalid PSA format")
            
            # T-Stage parsing
            t_raw = str(row["T-Stage"]).lower().strip()
            if pd.isna(row["T-Stage"]) or t_raw in ['nan', '']:
                error_reasons.append("Missing T-Stage")
                t_raw = None
            elif "tx" in t_raw:
                error_reasons.append("T-Stage unknown (Tx)")
                t_raw = None
            
            # N-Stage parsing
            n_raw = str(row["N-Stage"]).lower().strip()
            if pd.isna(row["N-Stage"]) or n_raw in ['nan', '']:
                error_reasons.append("Missing N-Stage")
                n_raw = None
            elif "nx" in n_raw:
                error_reasons.append("N-Stage unknown (Nx)")
                n_raw = None
            
            # M-Stage parsing
            m_raw = str(row.get("M-Stage", "")).lower().strip()
            if pd.isna(row.get("M-Stage")) or m_raw in ['nan', '']:
                m_raw = "m0"  # Assume M0 when M-Stage is missing
            elif m_raw == "not mentioned":
                m_raw = "m0"  # Assume M0 when noted as "not mentioned"
            elif "mx" in m_raw:
                error_reasons.append("M-Stage unknown (Mx)")
                m_raw = None
            
            # If any required information is missing, return Unknown
            if error_reasons:
                return f"Unknown ({', '.join(error_reasons)})"
            
            # AJCC classification logic
            if "m1" in m_raw:
                return "IVB"
            if "n1" in n_raw:
                return "IVA"
            if gg == 5:
                return "IIIC"
            if "t3" in t_raw or "t4" in t_raw:
                return "IIIB"
            if psa_val >= 20:
                return "IIIA"
            if gg == 1:
                if psa_val < 10 and ("t1" in t_raw or "t2a" in t_raw):
                    return "I"
                else:
                    return "IIA"
            if gg == 2:
                return "IIB"
            if gg in [3, 4]:
                return "IIC"
                
            return "Unknown (Classification logic error)"
            
        except Exception as e:
            return f"Unknown (Processing error: {str(e)})"

    # Consistency Validation
    
    def validate_consistency(self, row):
        """
        Comprehensive pathology staging consistency validation
        """
        inconsistencies = []
        
        # T-Stage vs Extraprostatic Extension
        if (row.get('T-Stage') in ['pT2a', 'pT2b', 'pT2c'] and 
            row.get('Extraprostatic Extension') == 'Present'):
            inconsistencies.append('T2_stage_with_EPE_should_be_T3a')
        
        # T-Stage vs Seminal Vesicle Invasion  
        if (row.get('T-Stage') in ['pT2a', 'pT2b', 'pT2c', 'pT3a'] and 
            row.get('Seminal Vesicle Invasion') == 'Present'):
            inconsistencies.append('SVI_present_should_be_T3b_or_higher')
        
        # N-Stage vs Lymph Node Metastasis Count
        self._validate_lymph_node_consistency(row, inconsistencies)
        
        # Gleason Score vs Grade Group
        self._validate_gleason_grade_consistency(row, inconsistencies)
        
        # Additional staging logic checks
        self._validate_additional_staging_logic(row, inconsistencies)
        
        return inconsistencies

    def _validate_lymph_node_consistency(self, row, inconsistencies):
        """Validate N-stage consistency with lymph node counts"""
        try:
            lymph_metastasis = row.get('Lymph Nodes with Metastasis')
            n_stage = row.get('N-Stage')
            
            if pd.notna(lymph_metastasis):
                lymph_count = float(str(lymph_metastasis).strip()) if str(lymph_metastasis).strip() else 0
                
                # N0 but positive nodes
                if n_stage == 'pN0' and lymph_count > 0:
                    inconsistencies.append('N0_stage_but_positive_lymph_nodes')
                
                # N1 but no positive nodes  
                elif n_stage == 'pN1' and lymph_count == 0:
                    inconsistencies.append('N1_stage_but_no_positive_lymph_nodes')
                    
        except (ValueError, TypeError):
            pass

    def _validate_gleason_grade_consistency(self, row, inconsistencies):
        """Validate Gleason score vs WHO Grade Group consistency"""
        try:
            primary = row.get('Primary Gleason Pattern')
            secondary = row.get('Secondary Gleason Pattern')
            grade = row.get('WHO Grade Group')
            
            if all(pd.notna(x) for x in [primary, secondary, grade]):
                p = int(float(str(primary).strip()))
                s = int(float(str(secondary).strip()))
                gleason_sum = p + s
                
                # WHO 2016 Grade Group mapping
                expected_grade = None
                if gleason_sum <= 6:
                    expected_grade = 'GG1'
                elif gleason_sum == 7 and p == 3:
                    expected_grade = 'GG2'
                elif gleason_sum == 7 and p == 4:
                    expected_grade = 'GG3'
                elif gleason_sum == 8:
                    expected_grade = 'GG4'
                elif gleason_sum >= 9:
                    expected_grade = 'GG5'
                
                if expected_grade and expected_grade != grade:
                    inconsistencies.append(f'Gleason_{gleason_sum}_inconsistent_with_{grade}')
                    
        except (ValueError, TypeError):
            pass

    def _validate_additional_staging_logic(self, row, inconsistencies):
        """Additional staging logic validations"""
        # T3b should have SVI
        if row.get('T-Stage') == 'pT3b' and row.get('Seminal Vesicle Invasion') == 'Absent':
            inconsistencies.append('T3b_stage_but_no_SVI')
        
        # T3a should have EPE
        if row.get('T-Stage') == 'pT3a' and row.get('Extraprostatic Extension') == 'Absent':
            inconsistencies.append('T3a_stage_but_no_EPE')

    # Knowledge Graph Construction
    
    def create_knowledge_graph(self, df: pd.DataFrame):
        """
        Build Neo4j knowledge graph with single-pass classification results
        
        Returns:
            list: Classification results for each patient including consistency validation
        """
        ajcc_results = []
        
        with self.driver.session() as session:
            for idx, row in df.iterrows():
                patient_id = str(row["ID"])
                
                # Perform single-pass classification
                ajcc_stage = self.calculate_ajcc_stage(row)
                
                # Perform consistency validation
                inconsistencies = self.validate_consistency(row)
                
                # Create graph nodes and relationships
                self._create_patient_nodes(session, patient_id, row, inconsistencies, ajcc_stage)
                
                # Store results for output
                result = {
                    "ID": row["ID"], 
                    "AJCC_Stage": ajcc_stage,
                    "Consistency_Errors": "; ".join(inconsistencies) if inconsistencies else "No errors",
                    "Error_Count": len(inconsistencies),
                    "Needs_Reextraction": len(inconsistencies) > 0
                }
                ajcc_results.append(result)
                
        return ajcc_results

    def _create_patient_nodes(self, session, patient_id, row, inconsistencies, ajcc_stage):
        """Create all nodes and relationships for a patient"""
        
        # Patient node
        session.run("""
            MERGE (p:Patient {id: $pid})
            SET p.original_id = $oid
        """, pid=patient_id, oid=row["ID"])

        # PathologyReport node
        report_id = f"report_{patient_id}"
        session.run("""
            MERGE (r:PathologyReport {id: $rid})
            SET r.histologic_subtype = $histologic_subtype,
                r.resection_margins = $resection_margins,
                r.lymph_nodes_examined = $lymph_nodes_examined,
                r.lymph_nodes_with_metastasis = $lymph_nodes_metastasis,
                r.primary_gleason = $primary_gleason,
                r.secondary_gleason = $secondary_gleason,
                r.tertiary_gleason = $tertiary_gleason,
                r.secondary_gleason_percentage = $secondary_percentage,
                r.extraprostatic_extension = $epe,
                r.seminal_vesicle_invasion = $svi,
                r.perineural_invasion = $pni,
                r.inconsistencies = $inconsistencies,
                r.needs_reextraction = $needs_reextraction
        """, 
        rid=report_id,
        histologic_subtype=row.get('Histologic Subtype', 'Not mentioned'),
        resection_margins=row.get('Resection Margins', 'Not mentioned'),
        lymph_nodes_examined=row.get('Number of Lymph Nodes examined'),
        lymph_nodes_metastasis=row.get('Lymph Nodes with Metastasis'),
        primary_gleason=row.get('Primary Gleason Pattern', 'Not mentioned'),
        secondary_gleason=row.get('Secondary Gleason Pattern', 'Not mentioned'),
        tertiary_gleason=row.get('Tertiary Gleason Pattern', 'Not mentioned'),
        secondary_percentage=row.get('Percentage of Secondary Gleason Pattern'),
        epe=row.get('Extraprostatic Extension', 'Not mentioned'),
        svi=row.get('Seminal Vesicle Invasion', 'Not mentioned'),
        pni=row.get('Perineural Invasion', 'Not mentioned'),
        inconsistencies=inconsistencies,
        needs_reextraction=len(inconsistencies) > 0)

        # Patient-Report relationship
        session.run("""
            MATCH (p:Patient {id: $pid}), (r:PathologyReport {id: $rid})
            MERGE (p)-[:HAS_REPORT]->(r)
        """, pid=patient_id, rid=report_id)

        # TNM Stage nodes
        self._create_tnm_nodes(session, report_id, row)
        
        # Grade Group node
        self._create_grade_node(session, report_id, row)
        
        # PSA Value node
        self._create_psa_node(session, report_id, row)
        
        # AJCC Stage node
        self._create_ajcc_node(session, report_id, ajcc_stage)

    def _create_tnm_nodes(self, session, report_id, row):
        """Create TNM stage nodes and relationships"""
        # M-Stage handling: assume M0 if "not mentioned"
        m_stage = row.get("M-Stage", "M0")
        if pd.isna(m_stage) or str(m_stage).strip().lower() == "not mentioned":
            m_stage = "M0"
        
        for stage_type, stage_value in [
            ("T", row["T-Stage"]),
            ("N", row["N-Stage"]),
            ("M", m_stage),
        ]:
            if not self._is_missing(stage_value):
                session.run("""
                    MERGE (s:TNM_Stage {code: $code, type: $stype})
                    SET s.value = $code
                """, code=stage_value, stype=stage_type)
                
                session.run(f"""
                    MATCH (r:PathologyReport {{id: $rid}})
                    MATCH (s:TNM_Stage {{code: $code, type: $stype}})
                    MERGE (r)-[:HAS_{stage_type}_STAGE]->(s)
                """, rid=report_id, code=stage_value, stype=stage_type)

    def _create_grade_node(self, session, report_id, row):
        """Create Grade Group node and relationship"""
        grade_group = row.get('WHO Grade Group')
        if not self._is_missing(grade_group):
            session.run("""
                MERGE (g:Grade {group: $grade_group})
            """, grade_group=grade_group)
            
            session.run("""
                MATCH (r:PathologyReport {id: $rid})
                MATCH (g:Grade {group: $grade_group})
                MERGE (r)-[:HAS_GRADE]->(g)
            """, rid=report_id, grade_group=grade_group)

    def _create_psa_node(self, session, report_id, row):
        """Create PSA Value node and relationship"""
        psa_value = row.get("Serum PSA")
        if not self._is_missing(psa_value):
            psa_id = f"psa_{psa_value}"
            session.run("""
                MERGE (p:PSA_Value {id: $pid})
                SET p.value = $val
            """, pid=psa_id, val=psa_value)
            
            session.run("""
                MATCH (r:PathologyReport {id: $rid})
                MATCH (p:PSA_Value {id: $pid})
                MERGE (r)-[:HAS_PSA]->(p)
            """, rid=report_id, pid=psa_id)

    def _create_ajcc_node(self, session, report_id, ajcc_stage):
        """Create AJCC Stage node and relationship"""
        session.run("MERGE (a:AJCC_Stage {stage: $s})", s=ajcc_stage)
        session.run("""
            MATCH (r:PathologyReport {id: $rid})
            MATCH (a:AJCC_Stage {stage: $s})
            MERGE (r)-[:CLASSIFIED_AS]->(a)
        """, rid=report_id, s=ajcc_stage)

    def create_ajcc_rules(self):
        """Add AJCC classification rules to the graph"""
        with self.driver.session() as session:
            rules = [
                # Stage I: T1-T2a AND N0 AND M0 AND GG1 AND PSA<10
                {"rule_id": "rule_stage_I", 
                 "condition": "T1-T2a AND N0 AND M0 AND Grade Group 1 AND PSA < 10", 
                 "result": "I"},

                # Stage IIA: T1-T2c AND N0 AND M0 AND ((GG1 AND PSA 10-19.9) OR other specified conditions)
                {"rule_id": "rule_stage_IIA", 
                 "condition": "T1-T2c AND N0 AND M0 AND Grade Group 1 AND PSA >= 10 (but < 20)", 
                 "result": "IIA"},

                # Stage IIB: T1-T2c AND N0 AND M0 AND GG2
                {"rule_id": "rule_stage_IIB", 
                 "condition": "T1-T2c AND N0 AND M0 AND Grade Group 2", 
                 "result": "IIB"},

                # Stage IIC: T1-T2c AND N0 AND M0 AND GG3-4
                {"rule_id": "rule_stage_IIC", 
                 "condition": "T1-T2c AND N0 AND M0 AND Grade Group 3-4", 
                 "result": "IIC"},

                # Stage IIIA: T1-T2 AND N0 AND M0 AND PSA >= 20
                {"rule_id": "rule_stage_IIIA", 
                 "condition": "T1-T2 AND N0 AND M0 AND PSA >= 20", 
                 "result": "IIIA"},

                # Stage IIIB: T3a AND N0 AND M0
                {"rule_id": "rule_stage_IIIB", 
                 "condition": "T3a AND N0 AND M0", 
                 "result": "IIIB"},

                # Stage IIIC: (T3b-T4 AND N0 AND M0) OR (Any T AND N0 AND M0 AND GG5)
                {"rule_id": "rule_stage_IIIC", 
                 "condition": "T3b-T4 AND N0 AND M0 OR Grade Group 5 AND N0 AND M0", 
                 "result": "IIIC"},

                # Stage IVA: Any T AND N1 AND M0
                {"rule_id": "rule_stage_IVA", 
                 "condition": "Any T AND N1 AND M0", 
                 "result": "IVA"},

                # Stage IVB: Any T AND Any N AND M1
                {"rule_id": "rule_stage_IVB", 
                 "condition": "Any T AND Any N AND M1", 
                 "result": "IVB"},
            ]

            for rule in rules:
                session.run("""
                    MERGE (r:Rule {id: $rule_id})
                    SET r.condition = $condition, r.result = $result
                """, **rule)

    # Analytics and Reporting
    
    def query_patient_info(self, patient_id):
        """Query detailed information for a specific patient"""
        with self.driver.session() as session:
            result = session.run("""
                MATCH (p:Patient {id: $patient_id})-[:HAS_REPORT]->(r:PathologyReport)
                OPTIONAL MATCH (r)-[:HAS_T_STAGE]->(t:TNM_Stage)
                OPTIONAL MATCH (r)-[:HAS_N_STAGE]->(n:TNM_Stage)
                OPTIONAL MATCH (r)-[:HAS_M_STAGE]->(m:TNM_Stage)
                OPTIONAL MATCH (r)-[:HAS_PSA]->(psa:PSA_Value)
                OPTIONAL MATCH (r)-[:HAS_GRADE]->(g:Grade)
                OPTIONAL MATCH (r)-[:CLASSIFIED_AS]->(a:AJCC_Stage)
                RETURN p, r, t.code as t_stage, n.code as n_stage, m.code as m_stage, 
                       psa.value as psa_value, g.group as grade_group, a.stage as ajcc_stage,
                       r.inconsistencies as inconsistencies, r.needs_reextraction as needs_reextraction
            """, patient_id=patient_id)
            return result.single()

    def find_inconsistent_cases(self):
        """Find all cases with consistency errors"""
        with self.driver.session() as session:
            result = session.run("""
                MATCH (p:Patient)-[:HAS_REPORT]->(r:PathologyReport)
                WHERE r.needs_reextraction = true
                RETURN p.id as patient_id, r.inconsistencies as inconsistencies
            """)
            return [(record["patient_id"], record["inconsistencies"]) for record in result]

    def get_statistics(self):
        """Get comprehensive knowledge graph statistics"""
        with self.driver.session() as session:
            # Basic counts
            patient_count = session.run("MATCH (p:Patient) RETURN count(p) as count").single()["count"]
            inconsistent_count = session.run("""
                MATCH (r:PathologyReport) WHERE r.needs_reextraction = true
                RETURN count(r) as count
            """).single()["count"]
            
            # AJCC Stage distribution
            ajcc_distribution = session.run("""
                MATCH (a:AJCC_Stage)<-[:CLASSIFIED_AS]-(r:PathologyReport)
                RETURN a.stage as stage, count(r) as count
                ORDER BY a.stage
            """)
            
            # Top consistency errors
            error_distribution = session.run("""
                MATCH (r:PathologyReport) WHERE r.needs_reextraction = true
                UNWIND r.inconsistencies as error
                RETURN error, count(error) as count
                ORDER BY count DESC
            """)
            
            return {
                'total_patients': patient_count,
                'ajcc_distribution': [(record["stage"], record["count"]) for record in ajcc_distribution],
                'inconsistent_cases': inconsistent_count,
                'top_consistency_errors': [(record["error"], record["count"]) for record in error_distribution]
            }

    def export_consistency_errors(self, df_with_classification: pd.DataFrame, output_file_path: str):
        """
        Export detailed consistency error analysis to Excel with multiple sheets
        
        Args:
            df_with_classification: Main dataset with classification results
            output_file_path: Excel file path
        """
        
        # Create error analysis DataFrame
        error_cases = df_with_classification[df_with_classification['Needs_Reextraction'] == True].copy()
        
        # Detailed error breakdown
        error_details = []
        for idx, row in error_cases.iterrows():
            if row['Consistency_Errors'] != "No errors":
                errors = row['Consistency_Errors'].split('; ')
                for error in errors:
                    error_details.append({
                        'Patient_ID': row['ID'],
                        'Error_Type': error,
                        'T_Stage': row.get('T-Stage', 'N/A'),
                        'N_Stage': row.get('N-Stage', 'N/A'),
                        'M_Stage': row.get('M-Stage', 'N/A'),
                        'Grade_Group': row.get('WHO Grade Group', 'N/A'),
                        'PSA': row.get('Serum PSA', 'N/A'),
                        'AJCC_Stage': row.get('AJCC_Stage', 'N/A')
                    })
        
        error_details_df = pd.DataFrame(error_details)
        
        # Error summary statistics
        if len(error_details_df) > 0:
            error_summary = error_details_df['Error_Type'].value_counts().reset_index()
            error_summary.columns = ['Error_Type', 'Count']
            error_summary['Percentage'] = (error_summary['Count'] / len(error_cases) * 100).round(2)
        else:
            error_summary = pd.DataFrame({'Error_Type': ['No errors found'], 'Count': [0], 'Percentage': [0]})
        
        # Write to Excel with multiple sheets
        with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
            # Main results sheet
            df_with_classification.to_excel(writer, sheet_name='Classification_Results', index=False)
            
            # Error cases only
            error_cases.to_excel(writer, sheet_name='Cases_With_Errors', index=False)
            
            # Detailed error breakdown
            if len(error_details_df) > 0:
                error_details_df.to_excel(writer, sheet_name='Error_Details', index=False)
            
            # Error summary
            error_summary.to_excel(writer, sheet_name='Error_Summary', index=False)
            
            # Statistics summary
            stats = self.get_statistics()
            stats_df = pd.DataFrame([
                ['Total Patients', stats['total_patients']],
                ['Cases with Consistency Errors', stats['inconsistent_cases']],
                ['Error Rate (%)', round(stats['inconsistent_cases']/stats['total_patients']*100, 2) if stats['total_patients'] > 0 else 0]
            ], columns=['Metric', 'Value'])
            
            # Add AJCC distribution
            ajcc_stats = pd.DataFrame(stats['ajcc_distribution'], columns=['AJCC_Stage', 'Count'])
            
            stats_df.to_excel(writer, sheet_name='Statistics', index=False, startrow=0)
            ajcc_stats.to_excel(writer, sheet_name='Statistics', index=False, startrow=len(stats_df)+2, startcol=0)
            
        print(f"Detailed consistency error analysis saved to '{output_file_path}'")
        print(f"Total cases with errors: {len(error_cases)}")
        print(f"Error rate: {len(error_cases)/len(df_with_classification)*100:.1f}%")
        
        return error_cases, error_details_df, error_summary


# Function to add AJCC Stage column to dataframe

def add_ajcc_stage_column(df):
    """
    Function to add AJCC_Stage column to dataframe
    """
    kg = ProstateCancerKnowledgeGraph()
    
    # Create AJCC_Stage column
    df['AJCC_Stage'] = df.apply(kg.calculate_ajcc_stage, axis=1)
    
    return df


# Main Execution Functions...

def _display_analytics(kg, df_with_classification):
    """Display comprehensive analytics including consistency errors"""
    print("\n" + "="*60)
    print("KNOWLEDGE GRAPH ANALYTICS REPORT")
    print("="*60)
    
    stats = kg.get_statistics()
    
    # Basic Statistics
    print(f"\nBASIC STATISTICS:")
    print(f"   Total patients: {stats['total_patients']}")
    print(f"   Cases with consistency errors: {stats['inconsistent_cases']}")
    error_rate = stats['inconsistent_cases']/stats['total_patients']*100 if stats['total_patients'] > 0 else 0
    print(f"   Error rate: {error_rate:.1f}%")
    
    # AJCC Stage Distribution
    print(f"\nAJCC STAGE DISTRIBUTION:")
    if stats['ajcc_distribution']:
        for stage, count in stats['ajcc_distribution']:
            percentage = count / stats['total_patients'] * 100 if stats['total_patients'] > 0 else 0
            print(f"   {stage:8}: {count:4d} cases ({percentage:5.1f}%)")
    else:
        print("   No classification data available")
    
    # Classification Failure Analysis
    print(f"\nCLASSIFICATION FAILURE ANALYSIS:")
    unknown_cases = df_with_classification[
        df_with_classification['AJCC_Stage'].str.startswith('Unknown', na=False)
    ]
    if len(unknown_cases) > 0:
        print(f"   Total Unknown cases: {len(unknown_cases)}")
        unknown_reasons = unknown_cases['AJCC_Stage'].value_counts()
        for reason, count in unknown_reasons.head(5).items():
            percentage = count / len(df_with_classification) * 100
            print(f"   {reason}: {count} cases ({percentage:.1f}%)")
    else:
        print("   All cases successfully classified!")
    
    # Top Consistency Errors
    print(f"\nTOP CONSISTENCY ERRORS:")
    if stats['top_consistency_errors']:
        for i, (error, count) in enumerate(stats['top_consistency_errors'][:5], 1):
            percentage = count / stats['inconsistent_cases'] * 100 if stats['inconsistent_cases'] > 0 else 0
            print(f"   {i}. {error}: {count} cases ({percentage:.1f}%)")
    else:
        print("   No consistency errors found!")
    
    print("\n" + "="*60)


def _create_sample_data():
    """Create sample data for testing when file is not found"""
    return pd.DataFrame({
        'ID': [1, 2, 3, 4],
        'T-Stage': ['pT2a', '', 'pT2c', 'Not mentioned'],
        'N-Stage': ['pN0', 'pN0', 'pN1', 'pN0'],
        'M-Stage': ['M0', '', 'M0', 'M0'],
        'Serum PSA': ['7.5', '25.2', '15.0', '12.0'],
        'WHO Grade Group': ['GG1', 'GG3', 'GG2', ''],
        'Extraprostatic Extension': ['Absent', 'Present', 'Absent', 'Absent'],
        'Seminal Vesicle Invasion': ['Absent', 'Absent', 'Absent', 'Absent'],
        'Primary Gleason Pattern': [3, 4, 3, 3],
        'Secondary Gleason Pattern': [3, 3, 4, 4],
        'Lymph Nodes with Metastasis': [0, 0, 2, 0]
    })


def main(excel_file_path: str = "input_data.xlsx", 
         output_file_path: str = "ajcc_classification_result.xlsx",
         log_file_path: str = None):
    """
    Main execution function for Prostate Cancer AJCC Classification
    
    Args:
        excel_file_path: Input Excel file path
        output_file_path: Output Excel file path with classification results and consistency errors
        log_file_path: Log file path (auto-generated if None)
    """
    
    # Logging Setup
    if log_file_path is None:
        timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file_path = f"ajcc_processing_log_{timestamp}.txt"
    
    tee = setup_logging(log_file_path)
    
    try:
        print("PROSTATE CANCER AJCC CLASSIFICATION SYSTEM")
        print("=" * 60)
        print(f"Processing started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Input file: {excel_file_path}")
        print(f"Output file: {output_file_path}")
        print(f"Log file: {log_file_path}")
        print("=" * 60)
        
        # Data Loading
        print("\nSTEP 1: DATA LOADING")
        print("-" * 30)
        
        try:
            df = pd.read_excel(excel_file_path)
            print(f"Excel file '{excel_file_path}' loaded successfully!")
            print(f"Total records: {len(df)}")
            print(f"Total columns: {len(df.columns)}")
            print(f"Columns: {list(df.columns)}")
            
            # Data quality check
            print(f"\nDATA QUALITY CHECK:")
            missing_data = df.isnull().sum()
            critical_cols = ['ID', 'T-Stage', 'N-Stage', 'WHO Grade Group', 'Serum PSA']
            for col in critical_cols:
                if col in df.columns:
                    missing_count = missing_data[col]
                    if col != 'ID':
                        not_mentioned_count = (df[col].astype(str).str.strip().str.lower().isin(['not mentioned', '', '-', 'nan', 'null'])).sum()
                        total_missing = missing_count + not_mentioned_count
                        print(f"   {col:20}: {missing_count:4d} null + {not_mentioned_count:4d} not_mentioned = {total_missing:4d} total missing ({total_missing/len(df)*100:5.1f}%)")
                    else:
                        print(f"   {col:20}: {missing_count:4d} missing ({missing_count/len(df)*100:5.1f}%)")
                else:
                    print(f"   {col:20}: Column not found!")
            
            print("\nData preview (first 3 rows):")
            print(df.head(3).to_string())
            
        except FileNotFoundError:
            print(f"File '{excel_file_path}' not found. Using sample data...")
            df = _create_sample_data()
            print(f"Sample data created with {len(df)} records")
        except Exception as e:
            print(f"Error loading file: {e}")
            return

        # Knowledge Graph Construction
        print(f"\nSTEP 2: CLASSIFICATION & KNOWLEDGE GRAPH CONSTRUCTION")
        print("-" * 55)
        
        kg = ProstateCancerKnowledgeGraph()
        
        print("Initializing knowledge graph...")
        kg.clear_database()
        print("   Database cleared")
        kg.create_constraints()
        print("   Constraints created")
        
        print(f"\nBuilding knowledge graph with classification...")
        ajcc_results = kg.create_knowledge_graph(df)
        print("   Patient nodes and relationships created")
        
        kg.create_ajcc_rules()
        print("   AJCC classification rules added")
        print("Knowledge graph construction completed!")

        # Results Export
        print(f"\nSTEP 3: RESULTS EXPORT")
        print("-" * 30)
        
        print("Preparing comprehensive results...")
        ajcc_df = pd.DataFrame(ajcc_results)
        df_with_classification = df.merge(ajcc_df, on=['ID'], how='left')
        
        print("Exporting results with error analysis...")
        error_cases, error_details, error_summary = kg.export_consistency_errors(
            df_with_classification, output_file_path
        )
        
        # Analytics and Reporting
        print(f"\nSTEP 4: ANALYTICS & REPORTING")
        print("-" * 35)
        
        _display_analytics(kg, df_with_classification)
        
        # Detailed Error Case Analysis
        print(f"\nDETAILED ERROR CASE ANALYSIS")
        print("-" * 35)
        
        if len(error_cases) > 0:
            print(f"Found {len(error_cases)} cases with consistency errors")
            print(f"Error rate: {len(error_cases)/len(df_with_classification)*100:.1f}%")
            
            # Show sample error cases
            print(f"\nSAMPLE ERROR CASES (first 5):")
            sample_errors = error_cases[['ID', 'Consistency_Errors', 'AJCC_Stage']].head(5)
            for i, (_, row) in enumerate(sample_errors.iterrows(), 1):
                print(f"   {i}. Patient {row['ID']}:")
                print(f"      Stage: {row['AJCC_Stage']}")
                print(f"      Errors: {row['Consistency_Errors']}")
                print()
            
            # Error type breakdown
            if len(error_details) > 0:
                print(f"ERROR TYPE BREAKDOWN:")
                error_type_counts = error_details['Error_Type'].value_counts()
                for error_type, count in error_type_counts.head(10).items():
                    percentage = count / len(error_details) * 100
                    print(f"   {error_type:40}: {count:3d} ({percentage:5.1f}%)")
            
        else:
            print("No consistency errors found in the dataset!")
            print("All pathology data is consistent with medical staging logic.")
        
        # Classification Success Analysis
        print(f"\nCLASSIFICATION SUCCESS ANALYSIS")
        print("-" * 40)
        
        classification_success = df_with_classification['AJCC_Stage'].apply(
            lambda x: 'Success' if not str(x).startswith('Unknown') else 'Failed'
        ).value_counts()
        
        for status, count in classification_success.items():
            percentage = count / len(df_with_classification) * 100
            print(f"   {status:10}: {count:4d} cases ({percentage:5.1f}%)")
        
        # Show failed cases
        failed_cases = df_with_classification[
            df_with_classification['AJCC_Stage'].str.startswith('Unknown', na=False)
        ]
        
        if len(failed_cases) > 0:
            print(f"\nFAILED CLASSIFICATION CASES:")
            for _, row in failed_cases.head(5).iterrows():
                print(f"   Patient {row['ID']}: {row['AJCC_Stage']}")
        
        # Final Summary
        print(f"\nPROCESSING SUMMARY")
        print("=" * 60)
        print(f"Successfully processed {len(df_with_classification)} patients")
        print(f"Results saved to: {output_file_path}")
        print(f"Log saved to: {log_file_path}")
        print(f"Processing completed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 60)
        
        success_rate = classification_success.get('Success', 0) / len(df_with_classification) * 100
        print(f"Classification success rate: {success_rate:.1f}%")
        
        if len(error_cases) == 0:
            print("Data quality: No consistency errors")
        else:
            quality_score = (1 - len(error_cases)/len(df_with_classification)) * 100
            print(f"Data quality score: {quality_score:.1f}%")
        
        print("\nAJCC classification processing completed successfully!")
        
    except Exception as e:
        print(f"\nCRITICAL ERROR during processing:")
        print(f"   Error: {e}")
        print(f"   Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
    finally:
        try:
            kg.close()
            print(f"\nDatabase connection closed")
        except:
            pass
        
        # Close logging
        tee.close()
        sys.stdout = tee.terminal  # Restore original stdout
        print(f"\nComplete log saved to: {log_file_path}")

if __name__ == "__main__":
    # Usage example:
    main('input_data.xlsx', 
         'kg_classification_result.xlsx',
         'kg_processing_log.txt')


