using Kavr
using CairoMakie

function plot(αmax, βmax, rmin, rmax, n, a, steps, θs, θo, B, βfluid, count)
    rh = 1 + √(1-a^2)

    
    f = Figure(resolution = (1000, 800))
    ax0 = Axis(f[1,1], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    xlims!(ax0, -αmax, αmax)
    ylims!(ax0, -βmax, βmax)

    ax1 = Axis(f[1,2], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    xlims!(ax1, -αmax, αmax)
    ylims!(ax1, -βmax, βmax)

    ax2 = Axis(f[2,1], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    xlims!(ax2, -αmax, αmax)
    ylims!(ax2, -βmax, βmax)
    
    ax3 = Axis(f[2,2], xlabel="α", ylabel="β",aspect=1, xlabelsize=20, ylabelsize=20, xticklabelsize=20, yticklabelsize=20)
    xlims!(ax3, -αmax, αmax)
    ylims!(ax3, -βmax, βmax)

    αvals = LinRange(-αmax,αmax, 2steps)
    βvals = LinRange(-βmax, βmax, 2steps)
    rsvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
    rsvals2 = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
    sinvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
    cosvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
    geotypevals = [zeros(2steps) for _ in 1:2steps]
    rootvals = [zeros(2steps) for _ in 1:2steps]

    #νθ = θo > θs 

    νθtrue = cos(θs)< abs(cos(θo)) ? (θo>θs) ⊻ (n%2==1) : false

    νθfalse = cos(θs)< abs(cos(θo)) ? (θo>θs) ⊻ (n%2==1) : true
    for _i in 1:(2steps)
        α = αvals[_i]
        for _j in 1:(2steps)
            β = βvals[_j]
            rstemp, νrtemp, numroots = rs(α, β, θs, θo, a, true, n)
            evalpol = (rstemp >= rh && rstemp != Inf)
            sintemp, costemp =  evalpol ? calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθtrue) : (0., 0.)

            rsvals[_i][_j] = rstemp
            rootvals[_i][_j] = numroots
            sinvals[_i][_j] = sintemp * αmax / 20
            cosvals[_i][_j] = costemp * βmax / 20

            rstemp, νrtemp, numroots = rs(α, β, θs, θo, a, false, n)
            evalpol = (rstemp >= rh && rstemp != Inf)
            sintemp, costemp =  evalpol ? calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθfalse) : (0., 0.)

            rsvals2[_i][_j] += rstemp
            rootvals[_i][_j] = numroots
            sinvals[_i][_j] += sintemp * αmax / 20
            cosvals[_i][_j] += costemp * βmax / 20

            geotypevals[_i][_j] = (Gθ(α, β, a, θs, θo, false,n)[2] ? 1 : 0)

        end
    end
    tickx = reduce(vcat,transpose.(sinvals)) 
    ticky = reduce(vcat,transpose.(cosvals))

    arrows!(ax2, .-αvals, βvals, .-tickx, ticky, arrowhead = ' ')
    if abs(cos(θo)) > cos(θs)
        rsvals .+= rsvals2
        contour!(ax1, .-αvals, βvals, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/10:rmax)
    else 
        contour!(ax1, .-αvals, βvals, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/10:rmax)
        contour!(ax1, .-αvals, βvals, reduce(vcat,transpose.(rsvals2)), levels=rmin:(rmax - rmin)/10:rmax)
    end
    heatmap!(ax0, .-αvals, βvals, reduce(vcat,transpose.(rootvals)))
    heatmap!(ax3, .-αvals, βvals, reduce(vcat,transpose.(geotypevals)))


    display(f)
end

a = -1.0
αmax = 10
βmax = 10
steps = 40
n = 0
θs = (60)*π/180
back = true

βv = 0.71
θz = π/2
ϕz = 90*π/180
ι = π/3
B = [sin(ι)*cos(ϕz), sin(ι)*sin(ϕz), cos(ι)]
βfluid = [βv, θz, ϕz]
rmin = 1 + √(1-a^2) 
rmax = 10

count = 0
@time for i in range( 1, 179, length=179)
    θo = i*π/180
    if θo == θs || (180-i)*π/180 == θs
        continue
    end
    plot(αmax, βmax, rmin, rmax, 0, a,steps, θs, θo, B, βfluid, count)
    count +=1
end

