rule download_ncbi:
     input:
          f"{config['outdir']}/logs/cleannamecheck.txt"
     params:
          required=config['download_genbank']['options']
     output:
          f"{config['outdir']}/logs/completed.txt"
     log: f"{config['outdir']}/logs/download_ncbi.log"
     threads:
          config['threads']
     run:
          if 'yes' in params.required:
               for gen in [line.rstrip() for line in open('genus.txt', 'r')]:
                    shell(f"ncbi-genome-download -g {gen} bacteria --parallel {threads} -v > {log}")
          shell(f"touch '{config['outdir']}/logs/completed.txt'")


rule gather_gbff:
     input:
          f"{config['outdir']}/logs/completed.txt"
     output:
          f"{config['outdir']}/logs/newroot.txt"
     params:
          required=config['download_genbank']['options']
     run:
          if 'yes' in params.required:
               files = glob.glob("**/*gbff.gz", recursive=True)
               outfile = open(f"{config['outdir']}/logs/newroot.txt", 'w')
               path=os.getcwd()
               for file in files:
                    print(Path(file).name)
                    shutil.copy(file, path)
                    outfile.write(f"{Path(file).name}\n")
          else:
               shell(f"touch '{config['outdir']}/logs/newroot.txt'")


rule gunzip_zipfiles:
     input:
          f"{config['outdir']}/logs/newroot.txt"
     output:
          f"{config['outdir']}/logs/gunzip_complete.txt"
     params:
          required=config['download_genbank']['options']
     run:
          if 'yes' in params.required:
               shell("gunzip *gbff.gz"),
               shell(f"touch '{config['outdir']}/logs/gunzip_complete.txt'")
          else:
               shell(f"touch '{config['outdir']}/logs/gunzip_complete.txt'")
