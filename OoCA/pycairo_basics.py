#!/usr/bin/env python

import cairo

def main():
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,250,250)
    ctx = cairo.Context(surface)

    ctx.set_source_rgba(1,1,1,1)
    ctx.rectangle(0,0,250,250)
    ctx.fill()
    ctx.stroke()

    ctx.set_source_rgba(0, 0, 1, 1)
    ctx.set_line_width(4)
    ctx.move_to(25,25)
    ctx.line_to(75,75)
    ctx.stroke()

    ctx.set_source_rgba(1, 0.5, 0, 1)
    ctx.rectangle(100,100,50,75)
    ctx.fill()
    ctx.stroke()


    surface.write_to_png("pycairo_basics.png")


if __name__ == "__main__":
    main()