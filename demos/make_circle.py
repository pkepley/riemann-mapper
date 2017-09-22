from svgpathtools import *

top_half = Arc(start=-3, radius=2.0+2j, rotation=0, large_arc=1, sweep=1, end=1)
bot_half = Arc(start=1.0, radius=2.0+2j, rotation=0, large_arc=1, sweep=1, end=-3)
circle = Path(top_half, bot_half)
wsvg(circle, filename='../imgs/circle.svg')


top_half = Arc(start=-3, radius=1.0+2j, rotation=0, large_arc=1, sweep=1, end=3)
bot_half = Arc(start=3.0, radius=1.0+2j, rotation=0, large_arc=1, sweep=1, end=-3)
ellipse = Path(top_half, bot_half)
wsvg(ellipse, filename='../imgs/ellipse.svg')


square = Path(Line(start=(0+0j), end=(1+0j)), Line(start=(1+0j), end=(1+1j)),
			  Line(start=(1+1j), end=(0+1j)), Line(start=(0+1j), end=(0+0j)))
wsvg(square, filename='../imgs/square.svg')
