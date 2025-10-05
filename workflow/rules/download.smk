# include: "globals.smk"


rule get_genomes_raw:
    input:
        IN_GENOMES,
    output:
        f"{RESULTS}/.genomes_raw.tsv",
    run:
        # weird, input/output substitution only works inside f-string
        utils.sort_filter_genomes(f"{input}", f"{output}", ONLY_REFSEQ)


rule get_metadata_raw:
    input:
        rules.get_genomes_raw.output,
    output:
        f"{RESULTS}/.genomes_metadata_raw.tsv",
    priority: 1
    retries: 3
    cache: True
    shell:
        """
        sed '1d' {input} | perl -ape '$_ = $F[1] . "\\n"' |\
        \
        datasets summary genome accession \
            --inputfile /dev/stdin \
            --as-json-lines |\
        tr -d '\\t' |\
        dataformat tsv genome |\
        tr -d '\\r' >| {output}
        """


rule get_metadata:
    input:
        rules.get_metadata_raw.output,
    output:
        f"{RESULTS}/genomes_metadata.tsv",
    shell:
        """
workflow/scripts/genome_metadata.R {input} >| {output}
"""


def get_genomes_dir(wc, output):
    return str(Path(output[0]).parent)


rule download_genomes:
    input:
        rules.get_genomes_raw.output,
    output:
        genomes=f"{RESULTS}/genomes/genomes.tsv",
        not_found=f"{RESULTS}/genomes/not_found.tsv",
    threads: workflow.cores
    params:
        genomes_dir=get_genomes_dir,
    shell:
        """
workflow/scripts/hydrate.py {threads} {params} {input}
"""


def params_output_name(wc, output):
    """
    Used by taxallnomy_targz
    """
    return str(Path(output[0]).name)


#rule taxallnomy_targz:
#    output:
#        f"{RESULTS}/taxallnomy.tar.gz",
#    priority: 1
#    retries: 3
#    cache: True
#    params:
#        url="https://sourceforge.net/projects/taxallnomy/files/latest/download",
#        output_name=params_output_name,
#    shell:
#        """
#        aria2c --dir {RESULTS}\
#            --continue=true --split 12\
#            --max-connection-per-server=16\
#            --min-split-size=1M\
#            --out={params.output_name}\
#            --quiet\
#            {params.url}
#        """
rule taxallnomy_targz:
    output:
        f"{RESULTS}/taxallnomy.tar.gz",
    priority: 1
    retries: 3
    cache: True
    params:
        url="https://sourceforge.net/projects/taxallnomy/files/latest/download",
        output_name=params_output_name,
    shell:
        r"""
        ( aria2c --dir {RESULTS} \
                 --continue=true --split=12 \
                 --max-connection-per-server=16 \
                 --min-split-size=1M \
                 --out={params.output_name} \
                 --quiet \
                 {params.url} \
          ) || \
        wget -O {RESULTS}/{params.output_name} {params.url}

        # sanity check: ensure we didn’t get an HTML page
        file {RESULTS}/{params.output_name} | grep -qi 'gzip compressed data'
        """


rule taxallnomy_linname:
    input:
        rules.taxallnomy_targz.output,
    output:
        f"{RESULTS}/taxallnomy_lin_name.tsv",
    cache: True
    params:
        ori=f"{RESULTS}/taxallnomy_database/taxallnomy_lin_name.tab",
    shell:
        """
tar --directory={RESULTS} -vxf {input}
mv {params.ori} {output}
"""


rule join_genomes_taxallnomy:
    input:
        taxallnomy=rules.taxallnomy_linname.output,
        #genomes=rules.get_metadata.output,
        genomes=rules.prepare_metadata.output,
    output:
        f"{RESULTS}/genomes_ranks.tsv",
    cache: True
    shell:
        """
workflow/scripts/cross.R {input} >| {output}
"""




rule local_genomes_stage:
    """
    Read a 4-col table (genome, faa, gff, tax_id) and stage files under
    results/genomes/<genome>/<genome>.faa|gff. Also write results/genomes/genomes.tsv
    compatible with the rest of the pipeline.
    Triggered only when LOCAL_GENOMES=True.
    """
    input:
        IN_GENOMES
    output:
        genomes=f"{RESULTS}/genomes/genomes.tsv",
        not_found=f"{RESULTS}/genomes/not_found.tsv",
    run:
        import csv, os, shutil
        from pathlib import Path

        in_table = Path(str(input[0]))
        out_dir = Path(f"{RESULTS}/genomes")
        out_dir.mkdir(parents=True, exist_ok=True)

        genomes_tsv = Path(str(output.genomes))
        not_found_tsv = Path(str(output.not_found))

        ok_rows = []
        bad_rows = []

        with open(in_table, newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            required = {"genome", "faa", "gff", "tax_id"}
            missing = required - set(r.fieldnames or [])
            if missing:
                raise RuntimeError(
                    f"Local genomes mode requires header with columns: {sorted(required)}\n"
                    f"Missing: {sorted(missing)} in {in_table}"
                )
            for row in r:
                genome = row["genome"].strip()
                faa = Path(row["faa"]).expanduser()
                gff = Path(row["gff"]).expanduser()
                tax_id = (row.get("tax_id") or "-1").strip()
                if not genome:
                    continue
                # check inputs
                if not faa.is_file() or not gff.is_file():
                    bad_rows.append((genome, str(faa), str(gff)))
                    continue

                # stage into results/genomes/<genome>/
                gdir = out_dir / genome
                gdir.mkdir(parents=True, exist_ok=True)

                # link if possible, else copy
                def link_or_copy(src, dst):
                    try:
                        if dst.exists() or dst.is_symlink():
                            dst.unlink()
                        os.symlink(src, dst)
                    except OSError:
                        shutil.copy2(src, dst)

                link_or_copy(faa, gdir / f"{genome}.faa")
                link_or_copy(gff, gdir / f"{genome}.gff")

                ok_rows.append((genome, tax_id))

        # Write genomes.tsv in the usual 4-column form (id, genome, refseq, version)
        with open(genomes_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["id", "genome", "refseq", "version"])
            for i, (genome, _) in enumerate(ok_rows, start=1):
                w.writerow([i, genome, False, 1])

        # Not-found report
        with open(not_found_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["id", "genome", "refseq", "version"])
            for i, (genome, faa, gff) in enumerate(bad_rows, start=1):
                w.writerow([i, f"{genome} (missing inputs: {faa} | {gff})", False, 1])


rule get_metadata_local:
    """
    Build genomes_metadata.tsv from the same local table so we can skip NCBI lookups.
    Only used when LOCAL_GENOMES=True.
    """
    input:
        IN_GENOMES
    output:
        f"{RESULTS}/genomes_metadata.tsv"
    run:
        import csv
        from pathlib import Path

        in_table = Path(str(input[0]))
        out_file = Path(str(output[0]))
        rows = []
        with open(in_table, newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                genome = row["genome"].strip()
                tax_id = (row.get("tax_id") or "-1").strip()
                if genome:
                    rows.append((genome, tax_id))
        with open(out_file, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["genome", "tax_id"])
            w.writerows(rows)


# Canonical targets the rest of the pipeline expects:
# - results/genomes/genomes.tsv
# - results/genomes_metadata.tsv

use rule download_genomes as _dl_genomes
use rule get_metadata     as _dl_metadata

rule prepare_genomes:
    output:
        genomes=f"{RESULTS}/genomes/genomes.tsv",
    input:
        _dl_genomes.output.genomes if not LOCAL_GENOMES else rules.local_genomes_stage.output.genomes

rule prepare_metadata:
    output:
        f"{RESULTS}/genomes_metadata.tsv",
    input:
        _dl_metadata.output if not LOCAL_GENOMES else rules.get_metadata_local.output

