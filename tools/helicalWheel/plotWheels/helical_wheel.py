import argparse
from plotWheels.helical_wheel import helical_wheel

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Helical Wheel")
    parser.add_argument("--sequence",dest="sequence",type=str)
    parser.add_argument("--seqRange",dest="seqRange",type=int,default=1)
    parser.add_argument("--t_size",dest="t_size",type=int,default=32)
    parser.add_argument("--rotation",dest="rotation",type=int,default=90)
    parser.add_argument("--output",dest="output",type=argparse.FileType("wb"), default="_helicalwheel.png")#dest="output",default="_helicalwheel.png")
    #### circle colors
    parser.add_argument("--f_A",dest="f_A", default="#ffcc33")
    parser.add_argument("--f_C",dest="f_C",default="#b5b5b5")
    parser.add_argument("--f_D",dest="f_D",default="#db270f")
    parser.add_argument("--f_E",dest="f_E",default="#db270f")
    parser.add_argument("--f_F",dest="f_F",default="#ffcc33")
    parser.add_argument("--f_G",dest="f_G",default="#b5b5b5")
    parser.add_argument("--f_H",dest="f_H",default="#12d5fc")
    parser.add_argument("--f_I",dest="f_I",default="#ffcc33")
    parser.add_argument("--f_K",dest="f_K",default="#12d5fc")
    parser.add_argument("--f_L",dest="f_L",default="#ffcc33")
    parser.add_argument("--f_M",dest="f_M",default="#ffcc33")
    parser.add_argument("--f_N",dest="f_N",default="#b5b5b5")
    parser.add_argument("--f_P",dest="f_P",default="#ffcc33")
    parser.add_argument("--f_Q",dest="f_Q",default="#b5b5b5")
    parser.add_argument("--f_R",dest="f_R",default="#12d5fc")
    parser.add_argument("--f_S",dest="f_S",default="#b5b5b5")
    parser.add_argument("--f_T",dest="f_T",default="#b5b5b5")
    parser.add_argument("--f_V",dest="f_V",default="#ffcc33")
    parser.add_argument("--f_W",dest="f_W",default="#ffcc33")
    parser.add_argument("--f_Y",dest="f_Y",default="#b5b5b5")
    ### text colors
    parser.add_argument("--t_A",dest="t_A",default="k")
    parser.add_argument("--t_C",dest="t_C",default="k")
    parser.add_argument("--t_D",dest="t_D",default="w")
    parser.add_argument("--t_E",dest="t_E",default="w")
    parser.add_argument("--t_F",dest="t_F",default="k")
    parser.add_argument("--t_G",dest="t_G",default="k")
    parser.add_argument("--t_H",dest="t_H",default="k")
    parser.add_argument("--t_I",dest="t_I",default="k")
    parser.add_argument("--t_K",dest="t_K",default="k")
    parser.add_argument("--t_L",dest="t_L",default="k")
    parser.add_argument("--t_M",dest="t_M",default="k")
    parser.add_argument("--t_N",dest="t_N",default="k")
    parser.add_argument("--t_P",dest="t_P",default="k")
    parser.add_argument("--t_Q",dest="t_Q",default="k")
    parser.add_argument("--t_R",dest="t_R",default="k")
    parser.add_argument("--t_S",dest="t_S",default="k")
    parser.add_argument("--t_T",dest="t_T",default="k")
    parser.add_argument("--t_V",dest="t_V",default="k")
    parser.add_argument("--t_W",dest="t_W",default="k")
    parser.add_argument("--t_Y",dest="t_Y",default="k")

    args = parser.parse_args()

    
    #print(type(args.output))

    f_colors = [args.f_A,args.f_C,args.f_D,args.f_E,args.f_F,args.f_G,args.f_H,args.f_I,args.f_K,
                args.f_L,args.f_M,args.f_N,args.f_P,args.f_Q,args.f_R,args.f_S,args.f_T,args.f_V,
                args.f_W,args.f_Y]

    t_colors = [args.t_A,args.t_C,args.t_D,args.t_E,args.t_F,args.t_G,args.t_H,args.t_I,args.t_K,
                args.t_L,args.t_M,args.t_N,args.t_P,args.t_Q,args.t_R,args.t_S,args.t_T,args.t_V,
                args.t_W,args.t_Y]
    
    colors = [f_colors, t_colors]

    tmp_file = "./tmp.png"

    helical_wheel(sequence=args.sequence,
                  colorcoding=colors[0],
                  text_color=colors[1],
                  seqRange=args.seqRange,
                  t_size=args.t_size,
                  rot=args.rotation,
                  filename=tmp_file,
                 )
    
    with open("tmp.png", "rb") as f:
        for line in f:
            args.output.write(line)
