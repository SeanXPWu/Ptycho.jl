abstract type ReconstructionMode end
struct All <: ReconstructionMode end

struct Row <: ReconstructionMode
    row::Integer
end


function list(mode::Row)

end

function list(mode::All)

end
