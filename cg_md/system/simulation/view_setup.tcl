mol delrep 0 top

# Protein
mol representation VDW 0.9 16
mol selection {name BB SC1 SC2 SC3 SC4}
mol color Chain
mol material GlassBubble
mol addrep top

# Water
#mol representation Points
#mol selection {name W}
#mol color ColorID 8
#mol addrep top

# IONS
mol representation VDW 0.3 8
mol selection {resname NA+ CL- ION}
mol color Name
mol material Transparent
mol addrep top

display projection Orthographic
#display depthcue on
#display ambientocclusion on
#display shadows on
