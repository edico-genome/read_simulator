from common import run_process
import gzip


def create_fastas(_inst):
    """
    use bcftools to create 2 modified fastas that can be used for read simulation
    """
    
    # only modify the fastas if there are variants in the VCFs
    found_a_variant = 0
    zipped = ".gz" in inst.module_settings["truth_set_vcf"]
    open_safe = gzip.open if zipped else open
    with open_safe(inst.module_settings["truth_set_vcf"]) as stream_in:
        for line in stream_in:
            if line[0] == "#":
                continue
            found_a_variant = 1
            break

    # use the unmodified fasta to simulate reads if there are no variants
    if not found_a_variant:
        inst.module_settings["mod_fasta_1"] = inst.module_settings["fasta_file"]
        inst.module_settings["mod_fasta_2"] = inst.module_settings["fasta_file"]
    else:
        # else create two modified fastas
        cmd_template = "bcftools consensus -c {liftover} -H {haplotype} " + \
                       "-f {fasta_file} {truth_set_vcf}"

        for hap in [1, 2]:
            mod_fasta = os.path.join(
                inst.module_settings["workdir"], "mod_fasta_{}.fa".format(hap))
            liftover = os.path.join(
                inst.module_settings["workdir"], "liftover_{}.txt".format(hap))
            options = {
                "liftover": liftover, 
                "haplotype": hap, 
                "fasta_file": inst.module_settings["fasta_file"],
                "truth_set_vcf": inst.module_settings['truth_set_vcf']}
            cmd = cmd_template.format(**options)
            try:
                run_process(cmd, inst.logger, mod_fasta)
            except:
                inst.logger.info(
                    "BCFTools failed. Let's continue and hope for the best")

            # update module settings
            inst.module_settings["mod_fasta_{}".format(hap)] = mod_fasta

