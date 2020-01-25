
def compare_taxa(tax1, tax2):
    score = 0
    joined = []
    for a, b in zip(tax1, tax2):
        if a == b:
            score += 1
            joined.append(a)
        else:
            break
    return(score, joined)

def print_lngs(tax1, tax2):
    """
    Takes two list of numbers and draws a simple line diagram showing the divergence
    """
    # try guessing the line width
    try:
        import shutil
        s = shutil.get_terminal_size()
        linewidth = s.columns
    except:
        linewidth = 80
    
    
    def makestr(lst):
        return([str(i) for i in lst])
    
    tax1 = makestr(tax1)
    tax2 = makestr(tax2)
    score, joined = compare_taxa(tax1, tax2)
    # shorten
    joined = makestr(joined)
    tax1 = tax1[score:]
    tax2 = tax2[score:]
    
    # lables
    label1 = "Model Lng."
    label2 = "Shared    "
    label3 = "Bin Lng.  "

    root = "{} {}".format(label2, " ".join(joined))
    tax1r = " ".join(tax1)
    tax2r = " ".join(tax2)
    # check if a lineage or root is to long
    def shorten(lst):
        nw = [lst[0]]
        for e in lst[1:-3]:
            if e != "...":
                nw.append(e)
        nw.append("...")
        nw.append(lst[-1])
        
        return(nw)
    
    # check if we can shorten the branches
    while len(root) + len(tax1r) > linewidth or len(root) + len(tax2r) > linewidth:
        canremove = 0
        # shorten first the branches then the root
        if len(root) + len(tax1r) > linewidth and len(tax1) >= 3:
            # keep track of if there is room to shorten
            if len(tax1) <= 3 and tax1[1] == "...":
                canremove = canremove - 1
            else:
                tax1 = shorten(tax1)
        else:
            canremove = canremove - 1
                
        if len(root) + len(tax2r) > linewidth and len(tax2) >= 3:
            # keep track of if there is room to shorten
            if len(tax2) <= 3 and tax2[1] == "...":
                canremove = canremove - 1
            else:
                tax2 = shorten(tax2)
        else:
            canremove = canremove - 1
        
        tax1r = " ".join(tax1)
        tax2r = " ".join(tax2)
        if canremove <= -2:
            break
            
    # check if we need to shorten the root
     # check if we can shorten the branches
    while len(root) + len(tax1r) > linewidth or len(root) + len(tax2r) > linewidth:
        canremove = 0
        # shorten first the branches then the root
        if len(root) + len(tax1r) > linewidth and len(joined) >= 3:
            # keep track of if there is room to shorten
            if len(joined) <= 3 and joined[1] == "...":
                canremove = canremove - 1
            else:
                joined = shorten(joined)
        else:
            canremove = canremove - 1

              
        root = "{} {}".format(label2, " ".join(joined))
        if len(root) + len(tax2r) > linewidth and len(joined) >= 3:
            # keep track of if there is room to shorten
            if len(joined) <= 3 and joined[1] == "...":
                canremove = canremove - 1
            else:
                joined = shorten(joined)
        else:
            canremove = canremove - 1

        root = "{} {}".format(label2, " ".join(joined))
        
        if canremove <= -1:
            break
     
                
    
    space = [" "] * (len(root) - len(label1))
    space = "".join(space)
    longspace = "".join([" "] * (len(root)))
    
    #print dendrogram
    print("")
    print("Infered lineage compared to the model lineage:")
    if len(tax1r) > 0:
        print(label1, end ="")
        print(space, end = "  ")
        print(tax1r)
        print(longspace, end = "/\n")
    else:
        print(label1)
        print("")
    print(root)
    if len(tax2r) > 0:
        print(longspace, end = "\\\n")
        print(label3, end ="")
        print(space, end = "  ")
        print(tax2r)
    else:
        print("")
        print(label3)
    print("\n")
