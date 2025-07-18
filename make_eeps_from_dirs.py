import os
from glob import glob
from mesa_reader import MesaData
import argparse

TEMPLATE = "input.template"

def write_input_file(tracks_dir, eeps_dir, isos_dir, tracks, Zval) -> str:

    t = MesaData(tracks[0])
    initialY = t.total_mass_he4[0] / t.star_mass[0]
    initialZ = 1 - ((t.toal_mass_he4[0] + t.total_mass_h1[0]) / t.star_mass[0])
    ntracks = len(tracks)
    tracks_list = ""
    for tr in tracks_list:
        tracks_list = f"{tracks_list}{tr}\n"

    with open(TEMPLATE, "r") as file:
        template_text = file.read()
    
    template_text = template_text.replace("<<Y>>", initialY)
    template_text = template_text.replace("<<Z>>", initialZ)
    template_text = template_text.replace("<<TRACKS>>", tracks_dir)
    template_text = template_text.replace("<<EEPS>>", eeps_dir)
    template_text = template_text.replace("<<ISOCHRONES>>", isos_dir)
    
    template_text = template_text.replace("<<NTRACKS>>", str(ntracks))
    template_text = template_text.replace("<<TRACK_FILES>>", tracks_list[:-1])

    return template_text

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tracks_dir", type=str, action="store")
    parser.add_argument("eeps_dir", type=str, action="store")
    parser.add_argument("isochrones_dir", type=str, action="store")
    parser.add_argument("Z", type=float, action="store")
    parser.add_argument("-g", "--glob_pattern", type=str, action="store", default="*.data")
    args = parser.parse_args()

    tracks_dir = args.tracks_dir
    pattern = args.glob_pattern
    eeps_dir = args.eeps_dir
    isos_dir = args.isochrones_dir
    Zval = float(args.Z)
    
    tracks = glob(os.path.join(tracks_dir, pattern))

    with open("input.run", "w") as file:
        file.write(write_input_file(tracks_dir, eeps_dir, isos_dir, tracks, None))
    
    