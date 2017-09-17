from svgpathtools import svg2paths, disvg, path

paths, attrs = svg2paths('../imgs/indiana_map.svg')

disvg(paths, stroke_widths=[5])
