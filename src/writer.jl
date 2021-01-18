export writexyz, writexyztraj

"""
write xyz file
"""
function writexyz(fn::String, g::Geometry)
    write(fn, xyzstring(g))
end

function xyzstring(g::Geometry)
    atoms = Dict(1=>"H", 5=>"B", 6=>"C", 15=>"P")
    text = "$(g.natom)\n"
    text *= "MolecularDynamics.jl\n"
    for (index, typeindex) in enumerate(g.vtype)
        atomnumber = typeindex
        atom = atoms[atomnumber]
        x, y, z = g.vr[:, index]
        text *= @sprintf("%s % 2.5f % 2.5f % 2.5f\n", atom, x, y, z)
    end
    return text
end

"""
write xyz file trajectory
"""
function writexyztraj(fn::String, vg::Vector{Geometry})
    text = ""
    for g in vg
        text *= xyzstring(g)
    end
    write(fn, text)
end