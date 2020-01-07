import os
import re
import xlrd
import pandas as pd


class Get_Specie_Parameters():

    def __init__(self):
        def get_my_path(bioinfo):
            bio_drive_dep = {}
            bio_drive_dep['af'] = bioinfo + "/mycobacteria_brucella/mycobacterium/tbc/af2122/script_dependents"
            bio_drive_dep['h37'] = bioinfo + "/mycobacteria_brucella/mycobacterium/tbc/h37/script_dependents"
            bio_drive_dep['ab1'] = bioinfo + "/mycobacteria_brucella/brucella/abortus1/script_dependents"
            bio_drive_dep['ab3'] = bioinfo + "/mycobacteria_brucella/brucella/abortus3/script_dependents"
            bio_drive_dep['suis1'] = bioinfo + "/mycobacteria_brucella/brucella/suis1/script_dependents"
            bio_drive_dep['suis2'] = bioinfo + "/mycobacteria_brucella/brucella/suis2/script_dependents"
            bio_drive_dep['suis3'] = bioinfo + "/mycobacteria_brucella/brucella/suis3/script_dependents"
            bio_drive_dep['mel1'] = bioinfo + "/mycobacteria_brucella/brucella/melitensis-bv1/script_dependents"
            bio_drive_dep['mel1b'] = bioinfo + "/mycobacteria_brucella/brucella/melitensis-bv1b/script_dependents"
            bio_drive_dep['mel2'] = bioinfo + "/mycobacteria_brucella/brucella/melitensis-bv2/script_dependents"
            bio_drive_dep['mel3'] = bioinfo + "/mycobacteria_brucella/brucella/melitensis-bv3/script_dependents"
            bio_drive_dep['canis'] = bioinfo + "/mycobacteria_brucella/brucella/canis/script_dependents"
            bio_drive_dep['ceti1'] = bioinfo + "/mycobacteria_brucella/brucella/ceti1/script_dependents"
            bio_drive_dep['ceti2'] = bioinfo + "/mycobacteria_brucella/brucella/ceti2/script_dependents"
            bio_drive_dep['ovis'] = bioinfo + "/mycobacteria_brucella/brucella/ovis/script_dependents"
            bio_drive_dep['neo'] = bioinfo + "/mycobacteria_brucella/brucella/neotomae/script_dependents"
            bio_drive_dep['para'] = bioinfo + "/mycobacterium/avium_complex/vsnp/NC_002944/script_dependents"
            bio_drive_dep['typhimurium-atcc13311'] = bioinfo + "/bacterial_identification/salmonella/vsnp/typhimurium-atcc13311/script_dependents"
            bio_drive_dep['typhimurium-14028S'] = bioinfo + "/bacterial_identification/salmonella/vsnp/typhimurium-14028S/script_dependents"
            bio_drive_dep['typhimurium-LT2'] = bioinfo + "/bacterial_identification/salmonella/vsnp/typhimurium-LT2/script_dependents"
            bio_drive_dep['heidelberg-SL476'] = bioinfo + "/bacterial_identification/salmonella/vsnp/heidelberg-SL476/script_dependents"
            bio_drive_dep['strepequi'] = bioinfo + "/bacterial_identification/streptococcus/vsnp/script_dependents"
            bio_drive_dep['infantis-FSIS1502169'] = bioinfo + "/bacterial_identification/salmonella/vsnp/infantis-FSIS1502169/dependents"
            bio_drive_dep['dublin-ATCC39184'] = bioinfo + "/bacterial_identification/salmonella/vsnp/dublin-ATCC39184/dependents"
            bio_drive_dep['blockley'] = bioinfo + "/bacterial_identification/salmonella/vsnp/blockley-79-1229/dependents"
            bio_drive_dep['kentucky-SA20030505'] = bioinfo + "/bacterial_identification/salmonella/vsnp/kentucky-SA20030505/dependents"
            bio_drive_dep['newport-USDA-ARS-USMARC-1925'] = bioinfo + "/bacterial_identification/salmonella/vsnp/newport-USDA-ARS-USMARC-1925/dependents"
            bio_drive_dep['senftenberg-NCTC10080'] = bioinfo + "/bacterial_identification/salmonella/vsnp/senftenberg-NCTC10080/dependents"
            bio_drive_dep['enteritidis-Durban'] = bioinfo + "/bacterial_identification/salmonella/vsnp/enteritidis-Durban/dependents"
            bio_drive_dep['montevideo-507440-20'] = bioinfo + "/bacterial_identification/salmonella/vsnp/montevideo-507440-20/dependents"
            bio_drive_dep['te_atcc35865'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_atcc35865/script_dependents"
            bio_drive_dep['te_09-0932'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_09-0932/script_dependents"
            bio_drive_dep['te_89-0490'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_89-0490/script_dependents"
            bio_drive_dep['te_92-0972'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_92-0972/script_dependents"
            bio_drive_dep['te_98-0554'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_98-0554/script_dependents"
            bio_drive_dep['te_mce9'] = bioinfo + "/bacterial_identification/taylorella/vsnp/te_mce9/script_dependents"
            bio_drive_dep['belize'] = bioinfo + "/Analysis/results/newcastle/vSNP/belize/script_dependents"
            bio_drive_dep['gua'] = bioinfo + "/Analysis/results/newcastle/vSNP/guatemala/script_dependents"
            return (bio_drive_dep)

        self.upload_to = None
        real_path = os.path.dirname(os.path.realpath(__file__))
        print("real path command --> {}".format(real_path))
        real_path = real_path.split('/')
        root_path = '/'.join(real_path)
        self.dependents_dir = root_path + "/dependencies"
        

    def choose(self, species_selection):

        def get_tb_codes():
            if os.path.isfile("/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/genotyping_codes.xlsx"):
                tb_geno_codes = ("/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/genotyping_codes.xlsx")
            elif os.path.isfile("/project/mycobacteria_brucella/mycobacterium/genotyping_codes.xlsx"):
                tb_geno_codes = ("/project/mycobacteria_brucella/mycobacterium/genotyping_codes.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/genotyping_codes.xlsx"):
            #     tb_geno_codes = ("/Users/tstuber/Desktop/to_delete/genotyping_codes.xlsx")
            else:
                return None
            wb = xlrd.open_workbook(tb_geno_codes)
            ws = wb.sheet_by_index(0)
            genotype_codes = {}
            for rownum in range(ws.nrows):
                new_name = str(ws.row_values(rownum)[0])
                new_name = new_name.rstrip()
                new_name = re.sub('[\/() ]', '_', new_name)
                new_name = re.sub('#', 'num', new_name)
                new_name = re.sub('_-', '_', new_name)
                new_name = re.sub('-_', '_', new_name)
                new_name = re.sub('__+', '_', new_name)
                new_name = re.sub('_$', '', new_name)
                new_name = re.sub('-$', '', new_name)
                new_name = re.sub(',', '', new_name)
                try:
                    elite_test = ws.row_values(rownum)[1]
                except IndexError:
                    #print("except IndexError: when changing names")
                    elite_test = ""
                #print("newname %s" % new_name)
                try:
                    if new_name[-1] != "_":
                        new_name = new_name + "_"
                except IndexError:
                    pass
                genotype_codes.update({new_name: elite_test})
            return genotype_codes

        def get_para_codes():
            if os.path.isfile("/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx"):
                tb_geno_codes = ("/Volumes/bioinfo/project/mycobacteria_brucella/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx")
            elif os.path.isfile("/project/mycobacteria_brucella/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx"):
                tb_geno_codes = ("/project/mycobacteria_brucella/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/genotyping_codes.xlsx"):
            #     tb_geno_codes = ("/Users/tstuber/Desktop/to_delete/genotyping_codes.xlsx")
            else:
                return None
            wb = xlrd.open_workbook(tb_geno_codes)
            ws = wb.sheet_by_index(0)
            genotype_codes = {}
            for rownum in range(ws.nrows):
                new_name = str(ws.row_values(rownum)[0])
                new_name = new_name.rstrip()
                new_name = re.sub('[\/() ]', '_', new_name)
                new_name = re.sub('#', 'num', new_name)
                new_name = re.sub('_-', '_', new_name)
                new_name = re.sub('-_', '_', new_name)
                new_name = re.sub('__+', '_', new_name)
                new_name = re.sub('_$', '', new_name)
                new_name = re.sub('-$', '', new_name)
                new_name = re.sub(',', '', new_name)
                try:
                    elite_test = ws.row_values(rownum)[1]
                except IndexError:
                    #print("except IndexError: when changing names")
                    elite_test = ""
                #print("newname %s" % new_name)
                try:
                    if new_name[-1] != "_":
                        new_name = new_name + "_"
                except IndexError:
                    pass
                genotype_codes.update({new_name: elite_test})
            return genotype_codes

        def get_brucella_codes():
            if os.path.isfile("/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
                bruc_geno_codes = ("/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx")
            elif os.path.isfile("/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
                bruc_geno_codes = ("/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx"):
            #     bruc_geno_codes = ("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx")
            else:
                return None
            print("Pulling in the Brucella genotype codes...")
            wb_in = xlrd.open_workbook(bruc_geno_codes)
            sheet_in = wb_in.sheet_by_index(1)
            genotype_codes = {}
            for row_data in sheet_in.col(32):
                row_data = row_data.value
                row_data = re.sub("/", "_", row_data)
                row_data = re.sub("\.", "_", row_data)
                row_data = re.sub("\*", "_", row_data)
                row_data = re.sub("\?", "_", row_data)
                row_data = re.sub("\(", "_", row_data)
                row_data = re.sub("\)", "_", row_data)
                row_data = re.sub("\[", "_", row_data)
                row_data = re.sub("\]", "_", row_data)
                row_data = re.sub(" ", "_", row_data)
                row_data = re.sub("{", "_", row_data)
                row_data = re.sub("}", "_", row_data)
                row_data = re.sub("\'", "_", row_data)
                row_data = re.sub("-_", "_", row_data)
                row_data = re.sub("_-", "_", row_data)
                row_data = re.sub("--", "_", row_data)
                row_data = re.sub("_$", "", row_data)
                row_data = re.sub("-$", "", row_data)
                row_data = re.sub("\'", "", row_data)
                row_data = re.sub(",", "", row_data)
                row_data = str(row_data)
                genotype_codes[row_data] = "" #the empty value can be used for elites
            return genotype_codes

        def get_heidelberg_codes():
            if os.path.isfile("/Volumes/bioinfo/project/bacterial_identification/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx"):
                heidel_geno_codes = ("/Volumes/bioinfo/project/bacterial_identification/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx")
            elif os.path.isfile("/project/bacterial_identification/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx"):
                heidel_geno_codes = ("/project/bacterial_identification/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx"):
            #     bruc_geno_codes = ("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx")
            else:
                return None

            in_df = pd.read_excel(heidel_geno_codes, index_col="Filename", usecols=[0, 1])
            in_dict = in_df.to_dict('dict')
            genotype_codes = in_dict['Isolate Name']
            return genotype_codes

        def get_vndv_codes():
            if os.path.isfile("/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx"):
                vndv_geno_codes = ("/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx")
            elif os.path.isfile("/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx"):
                vndv_geno_codes = ("/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx"):
            #     bruc_geno_codes = ("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx")
            else:
                return None

            in_df = pd.read_excel(vndv_geno_codes, index_col="Filename", usecols=[0, 1])
            in_dict = in_df.to_dict('dict')
            genotype_codes = in_dict['Isolate Name']
            return genotype_codes

        if species_selection == "montevideo-507440-20":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/montevideo-507440-20/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP007530.fasta",
                "gbk_file": [script_dependents + "/CP007530.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/montevideo-507440-20/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "enteritidis-Durban":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/enteritidis-Durban/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP007507.fasta",
                "gbk_file": [script_dependents + "/CP007507.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/enteritidis-Durban/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "blockley":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/blockley-79-1229/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/blockley-79-1229.fasta",
                "gbk_file": [script_dependents + "/blockley-79-1229.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/blockley-79-1229/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "newport-USDA-ARS-USMARC-1925":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/newport-USDA-ARS-USMARC-1925/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP025232.fasta",
                "gbk_file": [script_dependents + "/CP025232.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/newport-USDA-ARS-USMARC-1925/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "senftenberg-NCTC10080":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/senftenberg-NCTC10080/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/LS483465.fasta",
                "gbk_file": [script_dependents + "/LS483465.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/senftenberg-NCTC10080/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "kentucky-SA20030505":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/kentucky-SA20030505/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP022500.fasta",
                "gbk_file": [script_dependents + "/CP022500.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/kentucky-SA20030505/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "dublin-ATCC39184":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/dublin-ATCC39184/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP019179.fasta",
                "gbk_file": [script_dependents + "/CP019179.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/dublin-ATCC39184/tables_and_trees",
                "script_dependents": script_dependents,
            }

        elif species_selection == "infantis-FSIS1502169":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/infantis-FSIS1502169/dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/CP016406.fasta",
                "gbk_file": [script_dependents + "/CP016406.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/infantis-FSIS1502169/tables_and_trees",
                "script_dependents": script_dependents,
            }
        elif species_selection == "typhimurium-atcc13311":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/typhimurium-atcc13311/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP009102.fasta",
                "gbk_file": [script_dependents + "/NZ_CP009102.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/typhimurium-atcc13311/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "typhimurium-14028S":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/typhimurium-14028S/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/NC_016856-NC_016855.fasta",
                "hqs": script_dependents + "/highqualitysnps.vcf",
                "gbk_file": [script_dependents + "/NC_016856.gbk", script_dependents + "/NC_016855.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/typhimurium-14028S/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "typhimurium-LT2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/typhimurium-LT2/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/salmonella/typhimurium-LT2/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/AE006468.fasta",
                "gbk_file": [script_dependents + "/AE006468.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/typhimurium-LT2/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "heidelberg-SL476":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/salmonella/heidelberg-SL476/script_dependents"
            genotype_codes = get_heidelberg_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/salmonella/heidelberg-SL476/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_011083.fasta",
                "gbk_file": [script_dependents + "/NC_011083.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 300,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/salmonella/vsnp/heidelberg-SL476/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "strepequi":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/streptococcus/vsnp/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/streptococcus/vsnp/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_017582.fasta",
                "gbk_file": [script_dependents + "/NC_017582.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/streptococcus/vsnp/script2",
               "script_dependents": script_dependents,
            }
        elif species_selection == "te_atcc35865":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_atcc35865/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_atcc35865/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_018108.fasta",
                "gbk_file": [script_dependents + "/NC_018108.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_atcc35865/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_09-0932":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_09-0932/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_09-0932/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021201.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021201.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_09-0932/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_89-0490":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_89-0490/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_89-0490/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021199.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021199.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_89-0490/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_92-0972":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_92-0972/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_92-0972/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021060.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021060.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_92-0972/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_98-0554":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_98-0554/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_98-0554/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021246.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021246.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_98-0554/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_mce9":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bacterial_identification/taylorella/te_mce9/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_mce9/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_014914.fasta",
                "gbk_file": [script_dependents + "/NC_014914.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bacterial_identification/taylorella/vsnp/te_mce9/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ab1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/abortus1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/abortus1/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_006932-NC_006933.fasta",
                "gbk_file": [script_dependents + "/NC_006932.gbk", script_dependents + "/NC_006933.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                #N_threshold changed to 150 from 350 on 11Jan2019 for SNP calls within the vaccine strain.
                "N_threshold": 150,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/abortus1/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ab3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/abortus3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/abortus3/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007682-NZ_CP007683.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007682.gbk", script_dependents + "/NZ_CP007683.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/abortus3/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "canis":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/canis/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/canis/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_010103-NC_010104.fasta",
                "gbk_file": [script_dependents + "/NC_010103.gbk", script_dependents + "/NC_010104.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/canis/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ceti1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ceti1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/ceti1/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/Bceti1Cudo.fasta",
                "gbk_file": None,
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/ceti1/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ceti2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ceti2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/ceti2/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_022905-NC_022906.fasta",
                "gbk_file": [script_dependents + "/NC_022905.gbk", script_dependents + "/NC_022906.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/ceti2/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv1/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_003317-NC_003318.fasta",
                "gbk_file": [script_dependents + "/NC_003317.gbk", script_dependents + "/NC_003318.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv1/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel1b":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1b/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv1b/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP018508-NZ_CP018509.fasta",
                "gbk_file": [script_dependents + "/NZ_CP018508.gbk", script_dependents + "/NZ_CP018509.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv1b/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv2/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_012441-NC_012442.fasta",
                "gbk_file": [script_dependents + "/NC_012441.gbk", script_dependents + "/NC_012442.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv2/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv3/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007760-NZ_CP007761.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007760.gbk", script_dependents + "/NZ_CP007761.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/melitensis-bv3/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis1/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_017251-NC_017250.fasta",
                "gbk_file": [script_dependents + "/NC_017251.gbk", script_dependents + "/NC_017250.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis1/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis2/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_010169-NC_010167.fasta",
                "gbk_file": [script_dependents + "/NC_010169.gbk", script_dependents + "/NC_010167.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/Defining_SNPs.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis2/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis3/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007719-NZ_CP007718.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007719.gbk", script_dependents + "/NZ_CP007718.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis3/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis4":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis4/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis4/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/B-REF-BS4-40.fasta",
                "gbk_file": None,
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/suis4/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ovis":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ovis/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/ovis/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_009505-NC_009504.fasta",
                "gbk_file": [script_dependents + "/NC_009505.gbk", script_dependents + "/NC_009504.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/ovis/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "neo":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/neotomae/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/brucella/neotomae/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/KN046827.fasta",
                "gbk_file": [script_dependents + "/KN046827.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/brucella/neotomae/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "af":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                #if no species type in dictionary default to vSNP repo dependencies
                script_dependents = str(self.dependents_dir) + "/mycobacterium/tbc/af2122/script_dependents"
            genotype_codes = get_tb_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/mycobacterium/tbc/af2122/script1",
                "spoligo_db": script_dependents + "/spoligotype_db.txt",
                "reference": script_dependents + "/NC_002945v4.fasta",
                "gbk_file": [script_dependents + "/NC_002945v4.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx", # previous excelinfile
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/mycobacterium/tbc/af2122/script2", #previous bioinfoVCF
                "script_dependents": script_dependents,
            }
        elif species_selection == "h37":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/mycobacterium/tbc/h37/script_dependents"
            genotype_codes = get_tb_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacteria_brucella/mycobacterium/tbc/h37/script1",
                "spoligo_db": script_dependents + "/spoligotype_db.txt",
                "reference": script_dependents + "/NC_000962.fasta",
                "gbk_file": [script_dependents + "/NC_000962.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacteria_brucella/mycobacterium/tbc/h37/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "para-NC_002944":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/mycobacterium/avium_complex/NC_002944/script_dependents"
            genotype_codes = get_para_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_002944.fasta",
                "gbk_file": [script_dependents + "/NC_002944.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "para-CP033688":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/mycobacterium/avium_complex/CP033688/script_dependents"
            genotype_codes = get_para_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/CP033688/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/CP033688.fasta",
                "gbk_file": [script_dependents + "/CP033688.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/CP033688/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "flu":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/virus/h7n3/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": None, #str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/H7N3-reference.fasta",
                "gbk_file": [script_dependents + "/H7N3-Seq1.gbk", script_dependents + "/H7N3-Seq2.gbk", script_dependents + "/H7N3-Seq3.gbk", script_dependents + "/H7N3-Seq4.gbk", script_dependents + "/H7N3-Seq5.gbk", script_dependents + "/H7N3-Seq6.gbk", script_dependents + "/H7N3-Seq7.gbk", script_dependents + "/H7N3-Seq8.gbk", ],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 300,
                # "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": None, #str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "newcastle":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/virus/newcastle/script_dependents"
            genotype_codes = get_vndv_codes()
            parameters = {
                "upload_to": None, #str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/18-016505-001-fusion-HN.fasta",
                "gbk_file": None,
                "species": species_selection,
                "qual_threshold": 40,
                "N_threshold": 40,
                # "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": None, #str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "belize":
            #### Special for finding MKillian
            if os.path.isdir("/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian"):
                script_dependents = "/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/belize/script_dependents"
                upload_to = None
            elif os.path.isdir("/project/diagnostic_virology_laboratory/MKillian"):
                script_dependents = "/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/belize/script_dependents"
                upload_to = "/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/belize/step2"
            else:
                real_path = os.path.dirname(os.path.realpath(__file__))
                print("real path command --> {}".format(real_path))
                real_path = real_path.split('/')
                root_path = '/'.join(real_path)
                self.dependents_dir = root_path + "/dependencies"
            #####
            genotype_codes = get_vndv_codes()
            parameters = {
                "upload_to": None, #str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/KF767466.fasta",
                "gbk_file": [script_dependents + "/KF767466.gbk"],
                "species": species_selection,
                "qual_threshold": 40,
                "N_threshold": 40,
                # "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": upload_to,
                "script_dependents": script_dependents,
            }
        elif species_selection == "gua":
            #### Special for finding MKillian
            if os.path.isdir("/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian"):
                script_dependents = "/Volumes/bioinfo/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/guatemala/script_dependents"
                upload_to = None
            elif os.path.isdir("/project/diagnostic_virology_laboratory/MKillian"):
                script_dependents = "/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/guatemala/script_dependents"
                upload_to = "/project/diagnostic_virology_laboratory/MKillian/Analysis/results/newcastle/vSNP/guatemala/step2"
            else:
                real_path = os.path.dirname(os.path.realpath(__file__))
                print("real path command --> {}".format(real_path))
                real_path = real_path.split('/')
                root_path = '/'.join(real_path)
                self.dependents_dir = root_path + "/dependencies"
            #####
            genotype_codes = get_vndv_codes()
            parameters = {
                "upload_to": None, #str(self.upload_to) + "/mycobacterium/avium_complex/vsnp/NC_002944/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/Gua1_1407_2018.fasta",
                "gbk_file": [script_dependents + "/Gua1_1407_2018.gbk"],
                "species": species_selection,
                "qual_threshold": 40,
                "N_threshold": 40,
                # "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": upload_to,
                "script_dependents": script_dependents,
            }
        else:
            genotype_codes = None
            parameters = {
                "upload_to": None,
                "spoligo_db": None,
                "reference": None,
                "gbk_file": None,
                "species": None,
                "qual_threshold": None,
                "N_threshold": None,
                "definingSNPs": None,
                "remove_from_analysis": None,
                "filter_file": None,
                "step2_upload": None,
                "script_dependents": None,
            }

        return (parameters, genotype_codes)
