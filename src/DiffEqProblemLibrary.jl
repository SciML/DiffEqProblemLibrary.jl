module DiffEqProblemLibrary

import ODEProblemLibrary
import DAEProblemLibrary
import DDEProblemLibrary
import SDEProblemLibrary
import JumpProblemLibrary
import BVProblemLibrary
import NonlinearProblemLibrary

using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    @compile_workload begin
        isdefined(@__MODULE__, :ODEProblemLibrary)
        isdefined(@__MODULE__, :DAEProblemLibrary)
        isdefined(@__MODULE__, :DDEProblemLibrary)
        isdefined(@__MODULE__, :SDEProblemLibrary)
        isdefined(@__MODULE__, :JumpProblemLibrary)
        isdefined(@__MODULE__, :BVProblemLibrary)
        isdefined(@__MODULE__, :NonlinearProblemLibrary)
    end
end

end # module
