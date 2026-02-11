#!/usr/bin/env python

import cairo

def main():
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,1200,1500)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1,1,1)
    ctx.rectangle(0,0,1200,1500)
    ctx.fill()
    ctx.stroke()

    ctx.set_source_rgb(0.5,0.5,0.5)
    ctx.rectangle(189,
                  85,
                  193,
                  15)
    ctx.fill()
    ctx.stroke()


    surface.write_to_png("pycairo_basics.png")


if __name__ == "__main__":
    main()