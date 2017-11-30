
###########################################################
class PirsTruthSam(ModuleBase):
    default_settings = {}
    expected_settings = []
    promised_outputs = []

    def run(self):
        # truth sam
        self.module_settings["truth_sam"] = \
            os.path.join(self.module_settings["outdir"], "truth.sam")

        # xform path
        xform_script_path = os.path.join(this_dir_path, "alt_contig", "xform_pirs_read_info.pl")

        # get settings from previous runs
        for key in ["contig", "contig-1", "contig-2", "alt-sam", "read_info",
                    "contigs_combo_type", "contig-1-from", "contig-1-to", "contig-2-from",
                    "contig-2-to"]:
            if key in self.module_settings:
                self.module_settings[key] = self.module_settings[key]

        bcf_modified_name = {}
        for i in ["1", "2"]:
            _contig = "contig-{}".format(i)
            _contig_name = self.module_settings["contig-{}".format(i)]
            _from = self.module_settings["contig-{}-from".format(i)]
            _to = self.module_settings["contig-{}-to".format(i)]
            bcf_modified_name[_contig] = "{}:{}-{}".format(_contig_name, _from, _to)

        options = {
            "xform": xform_script_path,
            "read_info": self.module_settings["read_info"],
            "alt-sam": self.module_settings["alt-sam"],
            "bcf-contig-1": bcf_modified_name["contig-1"],
            "bcf-contig-2": bcf_modified_name["contig-2"],
            "orig-contig-1": self.module_settings["contig-1"],
            "orig-contig-2": self.module_settings["contig-2"],
            "contig-1-from": self.module_settings["contig-1-from"],
            "contig-2-from": self.module_settings["contig-2-from"],
            "truth_sam": self.module_settings["truth_sam"]
        }

        cmd = "{xform} --read-info <(zcat {read_info}) --alt-sam {alt-sam}"
        cmd += " --contig1 {bcf-contig-1} --contig1-rename {orig-contig-1}"
        cmd += " --contig1-offset {contig-1-from}"
        cmd += " --contig2 {bcf-contig-2} --contig2-rename {orig-contig-2}"
        cmd += " --contig2-offset {contig-2-from}"
        cmd += " > {truth_sam}"

        cmd = cmd.format(**options)

        try:
            self.logger.info("Truth SAM cmd: {}".format(cmd))
            # for redirected streams subprocess requires the actual bash executable
            res = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
            self.logger.info("Response: {}".format(res))
            self.logger.info("Truth sam created: {}".format(self.module_settings["truth_sam"]))
            self.db_api.upload_to_db('sam_gold', self.module_settings["truth_sam"])
        except Exception as e:
            raise PipelineExc('Failed to create gold Sam. Error message: %s' % e)
