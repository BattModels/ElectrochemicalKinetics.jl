using Unitful
import Unitful: @u_str

#default: do nothing
nounits_T(T::Real) = T
nounits_V(V::Real) = V

#array dispatches, to avoid generating new arrays
nounits_T(T::AbstractArray{E}) where E<:Real = T
nounits_V(V::AbstractArray{E}) where E<:Real = V

#transform to default units, if a correct unit is passed
nounits_T(T::Unitful.Temperature) = Unitful.ustrip(Unitful.uconvert(u"K",T))
nounits_V(V::Unitful.Voltage) = Unitful.ustrip(Unitful.uconvert(u"V",T))

#if an array of units is passed, transform to unitless array
nounits_T(T::AbstractArray{E}) where E<:Unitful.Temperature = nounits_T.(T)
nounits_V(V::AbstractArray{E}) where E<:Unitful.Voltage = nounits_V.(V)