# This file is part of StochasticGene.jl
# transition_rate_stacked_elements.jl

"""
    struct StackedElement

Structure for stacked R transition matrix elements
fields:
- `a`: row index
- `b`: column index
- `index`: rate vector index
- `pm`: sign of elements
- `buffer_size`: size of the buffer
- `buffer_state`: current state of the buffer
"""
struct StackedElement
    a::Int
    b::Int
    index::Int
    pm::Int
    buffer_size::Int
    buffer_state::Int
end

"""
    set_elements_stacked_R!(elements, R, buffer_size, nu::Vector)

Set elements for stacked R transition matrix
"""
function set_elements_stacked_R!(elements, R::Int, buffer_size::Int, nu::Vector)
    # First handle buffer accumulation
    for b = 0:buffer_size-1
        # Accumulation from buffer to next state
        push!(elements, StackedElement(b+1, b+2, nu[1], 1, buffer_size, b+1))
        push!(elements, StackedElement(b+2, b+2, nu[1], -1, buffer_size, b+2))
    end
    
    # Buffer ejection transitions
    for b = 1:buffer_size
        # Ejection to R states
        for r = 1:R
            push!(elements, StackedElement(b, r, nu[2], 1, buffer_size, b))
            push!(elements, StackedElement(r, b, nu[2], -1, buffer_size, b))
        end
    end
    
    # Regular R transitions
    for r = 1:R-1
        push!(elements, StackedElement(r, r+1, nu[3], 1, buffer_size, 0))
        push!(elements, StackedElement(r+1, r+1, nu[3], -1, buffer_size, 0))
    end
end

"""
    set_elements_stacked_R(R, buffer_size, nu::Vector)

Return stacked R transition matrix elements
"""
function set_elements_stacked_R(R::Int, buffer_size::Int, nu::Vector)
    elements = Vector{StackedElement}(undef, 0)
    set_elements_stacked_R!(elements, R, buffer_size, nu)
    elements
end
