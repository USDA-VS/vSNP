import os
import re
import xlrd
import pandas as pd


class Get_Specie_Parameters():

    def __init__(self):
        def get_my_path(bioinfo):
            bio_drive_dep = {}
            bio_drive_dep['af'] = bioinfo + "/mycobacterium/tbc/af2122/script_dependents"
            bio_drive_dep['h37'] = bioinfo + "/mycobacterium/tbc/h37/script_dependents"
            bio_drive_dep['ab1'] = bioinfo + "/brucella/abortus1/script_dependents"
            bio_drive_dep['ab3'] = bioinfo + "/brucella/abortus3/script_dependents"
            bio_drive_dep['suis1'] = bioinfo + "/brucella/suis1/script_dependents"
            bio_drive_dep['suis2'] = bioinfo + "/brucella/suis2/script_dependents"
            bio_drive_dep['suis3'] = bioinfo + "/brucella/suis3/script_dependents"
            bio_drive_dep['mel1'] = bioinfo + "/brucella/melitensis-bv1/script_dependents"
            bio_drive_dep['mel1b'] = bioinfo + "/brucella/melitensis-bv1b/script_dependents"
            bio_drive_dep['mel2'] = bioinfo + "/brucella/melitensis-bv2/script_dependents"
            bio_drive_dep['mel3'] = bioinfo + "/brucella/melitensis-bv3/script_dependents"
            bio_drive_dep['canis'] = bioinfo + "/brucella/canis/script_dependents"
            bio_drive_dep['ceti1'] = bioinfo + "/brucella/ceti1/script_dependents"
            bio_drive_dep['ceti2'] = bioinfo + "/brucella/ceti2/script_dependents"
            bio_drive_dep['ovis'] = bioinfo + "/brucella/ovis/script_dependents"
            bio_drive_dep['neo'] = bioinfo + "/brucella/neotomae/script_dependents"
            bio_drive_dep['para'] = bioinfo + "/mycobacterium/avium_complex/vsnp/NC_002944/script_dependents"
            bio_drive_dep['typhimurium-atcc13311'] = bioinfo + "/bi/salmonella/vsnp/typhimurium-atcc13311/script_dependents"
            bio_drive_dep['typhimurium-14028S'] = bioinfo + "/bi/salmonella/vsnp/typhimurium-14028S/script_dependents"
            bio_drive_dep['typhimurium-LT2'] = bioinfo + "/bi/salmonella/vsnp/typhimurium-LT2/script_dependents"
            bio_drive_dep['heidelberg-SL476'] = bioinfo + "/bi/salmonella/vsnp/heidelberg-SL476/script_dependents"
            bio_drive_dep['te_atcc35865'] = bioinfo + "/bi/taylorella/vsnp/te_atcc35865/script_dependents"
            bio_drive_dep['te_09-0932'] = bioinfo + "/bi/taylorella/vsnp/te_09-0932/script_dependents"
            bio_drive_dep['te_89-0490'] = bioinfo + "/bi/taylorella/vsnp/te_89-0490/script_dependents"
            bio_drive_dep['te_92-0972'] = bioinfo + "/bi/taylorella/vsnp/te_92-0972/script_dependents"
            bio_drive_dep['te_98-0554'] = bioinfo + "/bi/taylorella/vsnp/te_98-0554/script_dependents"
            bio_drive_dep['te_mce9'] = bioinfo + "/bi/taylorella/vsnp/te_mce9/script_dependents"
            return (bio_drive_dep)
        
        if os.path.isdir("/Volumes/root/TStuber/Results"):
            bioinfo = "/Volumes/root/TStuber/Results"
            self.bio_drive_dep = get_my_path(bioinfo)
            self.upload_to = None
        elif os.path.isdir("/bioinfo11/TStuber/Results"):
            bioinfo = "/bioinfo11/TStuber/Results/"
            self.bio_drive_dep = get_my_path(bioinfo)
            self.upload_to = "/bioinfo11/TStuber/Results"
        else:
            self.upload_to = None
        real_path = os.path.dirname(os.path.realpath(__file__))
        print("real path command --> {}".format(real_path))
        real_path = real_path.split('/')
        root_path = '/'.join(real_path)
        self.dependents_dir = root_path + "/dependencies"
        

    def choose(self, species_selection):

        def get_tb_codes():
            if os.path.isfile("/Volumes/root/TStuber/Results/mycobacterium/genotyping_codes.xlsx"):
                tb_geno_codes = ("/Volumes/root/TStuber/Results/mycobacterium/genotyping_codes.xlsx")
            elif os.path.isfile("/bioinfo11/TStuber/Results/mycobacterium/genotyping_codes.xlsx"):
                tb_geno_codes = ("/bioinfo11/TStuber/Results/mycobacterium/genotyping_codes.xlsx")
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
            if os.path.isfile("/Volumes/root/TStuber/Results/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx"):
                tb_geno_codes = ("/Volumes/root/TStuber/Results/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx")
            elif os.path.isfile("/bioinfo11/TStuber/Results/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx"):
                tb_geno_codes = ("/bioinfo11/TStuber/Results/mycobacterium/avium_complex/metadata/avium_genotyping_codes.xlsx")
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
            if os.path.isfile("/Volumes/root/TStuber/Results/bi/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx"):
                heidel_geno_codes = ("/Volumes/root/TStuber/Results/bi/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx")
            elif os.path.isfile("/bioinfo11/TStuber/Results/bi/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx"):
                heidel_geno_codes = ("/bioinfo11/TStuber/Results/bi/salmonella/metadata/Heidelberg_vSNP_Metatada.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx"):
            #     bruc_geno_codes = ("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx")
            else:
                return None

            in_df = pd.read_excel(heidel_geno_codes, index_col="Filename", usecols=[0, 1])
            in_dict = in_df.to_dict('dict')
            genotype_codes = in_dict['Isolate Name']
            return genotype_codes

        def get_vndv_codes():
            if os.path.isfile("/Volumes/root/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx"):
                vndv_geno_codes = ("/Volumes/root/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx")
            elif os.path.isfile("/bioinfo14/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx"):
                vndv_geno_codes = ("/bioinfo14/MKillian/Analysis/results/newcastle/vSNP/vsnp_metadata.xlsx")
            # elif os.path.isfile("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx"):
            #     bruc_geno_codes = ("/Users/tstuber/Desktop/to_delete/ALL_WGS.xlsx")
            else:
                return None

            in_df = pd.read_excel(vndv_geno_codes, index_col="Filename", usecols=[0, 1])
            in_dict = in_df.to_dict('dict')
            genotype_codes = in_dict['Isolate Name']
            return genotype_codes

        if species_selection == "typhimurium-atcc13311":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/salmonella/typhimurium-atcc13311/script_dependents"
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
                "step2_upload": str(self.upload_to) + "/bi/salmonella/vsnp/typhimurium-atcc13311/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "typhimurium-14028S":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/salmonella/typhimurium-14028S/script_dependents"
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
                "step2_upload": str(self.upload_to) + "/bi/salmonella/vsnp/typhimurium-14028S/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "typhimurium-LT2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/salmonella/typhimurium-LT2/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/salmonella/typhimurium-LT2/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/AE006468.fasta",
                "gbk_file": [script_dependents + "/AE006468.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/salmonella/vsnp/typhimurium-LT2/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "heidelberg-SL476":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/salmonella/heidelberg-SL476/script_dependents"
            genotype_codes = get_heidelberg_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/salmonella/heidelberg-SL476/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_011083.fasta",
                "gbk_file": [script_dependents + "/NC_011083.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 300,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/salmonella/vsnp/heidelberg-SL476/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_atcc35865":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_atcc35865/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_atcc35865/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_018108.fasta",
                "gbk_file": [script_dependents + "/NC_018108.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_atcc35865/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_09-0932":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_09-0932/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_09-0932/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021201.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021201.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_09-0932/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_89-0490":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_89-0490/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_89-0490/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021199.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021199.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_89-0490/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_92-0972":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_92-0972/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_92-0972/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021060.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021060.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_92-0972/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_98-0554":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_98-0554/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_98-0554/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP021246.fasta",
                "gbk_file": [script_dependents + "/NZ_CP021246.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_98-0554/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "te_mce9":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/bi/taylorella/te_mce9/script_dependents"
            genotype_codes = None
            parameters = {
                "upload_to": str(self.upload_to) + "/bi/taylorella/vsnp/te_mce9/script1",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_014914.fasta",
                "gbk_file": [script_dependents + "/NC_014914.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/bi/taylorella/vsnp/te_mce9/script2",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ab1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/abortus1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/abortus1/data",
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
                "step2_upload": str(self.upload_to) + "/brucella/abortus1/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ab3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/abortus3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/abortus3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007682-NZ_CP007683.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007682.gbk", script_dependents + "/NZ_CP007683.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/abortus3/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "canis":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/canis/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/canis/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_010103-NC_010104.fasta",
                "gbk_file": [script_dependents + "/NC_010103.gbk", script_dependents + "/NC_010104.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/canis/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ceti1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ceti1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ceti1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/Bceti1Cudo.fasta",
                "gbk_file": None,
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/ceti1/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ceti2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ceti2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ceti2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_022905-NC_022906.fasta",
                "gbk_file": [script_dependents + "/NC_022905.gbk", script_dependents + "/NC_022906.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/ceti2/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_003317-NC_003318.fasta",
                "gbk_file": [script_dependents + "/NC_003317.gbk", script_dependents + "/NC_003318.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv1/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel1b":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1b/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv1b/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP018508-NZ_CP018509.fasta",
                "gbk_file": [script_dependents + "/NZ_CP018508.gbk", script_dependents + "/NZ_CP018509.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv1b/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_012441-NC_012442.fasta",
                "gbk_file": [script_dependents + "/NC_012441.gbk", script_dependents + "/NC_012442.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv2/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "mel3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007760-NZ_CP007761.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007760.gbk", script_dependents + "/NZ_CP007761.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv3/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis1":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis1/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_017251-NC_017250.fasta",
                "gbk_file": [script_dependents + "/NC_017251.gbk", script_dependents + "/NC_017250.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/suis1/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis2":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis2/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_010169-NC_010167.fasta",
                "gbk_file": [script_dependents + "/NC_010169.gbk", script_dependents + "/NC_010167.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/Defining_SNPs.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/suis2/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis3":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis3/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007719-NZ_CP007718.fasta",
                "gbk_file": [script_dependents + "/NZ_CP007719.gbk", script_dependents + "/NZ_CP007718.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/suis3/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "suis4":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/suis4/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis4/data",
                "spoligo_db": None,
                "reference": script_dependents + "/B-REF-BS4-40.fasta",
                "gbk_file": None,
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/suis4/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "ovis":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/ovis/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ovis/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_009505-NC_009504.fasta",
                "gbk_file": [script_dependents + "/NC_009505.gbk", script_dependents + "/NC_009504.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/ovis/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "neo":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/brucella/neotomae/script_dependents"
            genotype_codes = get_brucella_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/neotomae/data",
                "spoligo_db": None,
                "reference": script_dependents + "/KN046827.fasta",
                "gbk_file": [script_dependents + "/KN046827.gbk"],
                "species": species_selection,
                "qual_threshold": 300,
                "N_threshold": 350,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx",
                "step2_upload": str(self.upload_to) + "/brucella/neotomae/vcfs",
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
                "upload_to": str(self.upload_to) + "/mycobacterium/tbc/af2122/script1",
                "spoligo_db": script_dependents + "/spoligotype_db.txt",
                "reference": script_dependents + "/NC_002945v4.fasta",
                "gbk_file": [script_dependents + "/NC_002945v4.gbk"],
                "species": species_selection,
                "qual_threshold": 150,
                "N_threshold": 150,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx", # previous excelinfile
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/af2122/script2", #previous bioinfoVCF
                "script_dependents": script_dependents,
            }
        elif species_selection == "h37":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/mycobacterium/tbc/h37/script_dependents"
            genotype_codes = get_tb_codes()
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/tbc/h37/script1",
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
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/h37/script2",
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
                "step2_upload": None, #str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/vcfs",
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
                "step2_upload": None, #str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/vcfs",
                "script_dependents": script_dependents,
            }
        elif species_selection == "belize":
            try:
                script_dependents = self.bio_drive_dep[species_selection]
            except (KeyError, AttributeError):
                script_dependents = str(self.dependents_dir) + "/virus/belize/script_dependents"
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
                "step2_upload": None, #str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/vcfs",
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
