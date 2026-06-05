import os

output = "ATP5F1A_combined_sequence.txt"
script_name = os.path.basename(__file__)   

with open(output, "w", encoding="utf-8") as outfile:
    for item in os.listdir("."):
    
        if (os.path.isfile(item) and
            item != script_name and
            item != output):
            with open(item, "r", encoding="utf-8", errors="ignore") as infile:
                outfile.write(infile.read())
            outfile.write("\n\n")