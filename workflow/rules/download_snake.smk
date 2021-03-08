rule download_ncbi:
      params:
           required=config['download_genbank']['options']
      output:
           "logs/completed.txt"
      threads:
           config['threads']
      run:
           if 'yes' in params.required:
               for gen in [line.rstrip() for line in open('genus.txt', 'r')]:
                   shell(f"ncbi-genome-download --genera {gen} bacteria --parallel {threads}")
           shell("touch 'logs/completed.txt'")


rule gather_gbff:
      input:
           "logs/completed.txt"
      output:
           "logs/newroot.txt"
      params:
           required=config['download_genbank']['options']
      run:
           if 'yes' in params.required:
               files = glob.glob("**/*gbff.gz", recursive=True)
	       outfile = open("logs/newroot.txt", 'w')
	       path=os.getcwd()
	       for file in files:
                   print(Path(file).name)
                   shutil.copy(file, path)
                   outfile.write(Path(file).name+"\n")
           else:
               shell("touch 'logs/newroot.txt'")


rule gunzip_zipfiles:
       input:
            "logs/newroot.txt"
       output:
            "logs/gunzip_complete.txt"
       params:
            required=config['download_genbank']['options']
       run:
            if 'yes' in params.required:
                shell("gunzip *gbff.gz")
            shell("touch 'logs/gunzip_complete.txt'")
