# Author : Damien Dufour, University of Oslo 

# config file 
configfile: "config/config.yaml"  # replace 

# Directories

MAIN_DIR_ABS = os.getcwd()
OUT_DIR  = os.path.join(config["results_dir"], config["title"])

try:
    os.mkdir(OUT_DIR)
except FileExistsError:
    print("Folder already exists!")


rule main:
    input:
        data = config["bw_file"]
    log:
        os.path.join(OUT_DIR, "logs", "track.log")
    params:
        title              = config["title"],\
        coordinates        = config["coordinates"], \
        max                = config["max"], \
        day                = config["day"], \
        modality           = config["modality"], \
        ttt                = config["ttt"], \
        genome             = config["genome"], \
        interaction        = config["interaction"], \
        interaction_file   = config["interaction_file"], \
        output_dir         = OUT_DIR
    output:
        data = os.path.join(OUT_DIR,"data/",config["title"]+"_input_data.csv"),\
        figures = os.path.join(OUT_DIR,"figures/",config["title"]+"_track.pdf")
    shell:
        """Rscript script/track_visualisation.R \
        --input {input.data} \
        --output_dir {params.output_dir} \
        --title {params.title}\
        --coordinates {params.coordinates}\
        --max {params.max}\
        --day {params.day}\
        --modality {params.modality}\
        --ttt {params.ttt}\
        --genome {params.genome}\
        --interaction {params.interaction}\
        --interaction_file {params.interaction_file} \
        > {log} 2>&1"""

