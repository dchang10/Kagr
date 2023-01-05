using Revise
using Kagr
using CairoMakie

prof = profile(bam, r)*redshift^(3+bam.α)

DblPower(r, r0, a, b) = (r/r0)^a/(1 + (r/r0)^(a+b))
function Kagr.calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθfalse, rpeak, p1, p2) 
    return Kagr.calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθfalse) .* DblPower(rstemp, rpeak, p1, p2)
end

function plot(αmax, βmax, rmin, rmax, ntemp, a, steps, θs, θo, B, βfluid, rpeak, p1, p2, pa)
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

    for n in 0:ntemp

        αvalstemp = LinRange(-αmax,αmax, 2steps)
        βvalstemp = LinRange(-βmax, βmax, 2steps)
        αvals = (cos(pa) .* αvalstemp) .+ (sin(pa) .* βvalstemp)
        βvals = (cos(pa) .* βvalstemp) .- (sin(pa) .* αvalstemp)

        rsvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        rsvals2 = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        rsvals3 = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        rsvals4 = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        sinvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        cosvals = [zeros(2steps) for _ in 1:2steps]#[Vector{Float64}(0.0 ,2steps) for _ in 1:2steps]
        geotypevals = [zeros(2steps) for _ in 1:2steps]
        rootvals = [zeros(2steps) for _ in 1:2steps]

        #νθ = θo > θs 

        isincone = abs(cos(θs)) < abs(cos(θo))
        νθtrue =  isincone ? (n%2==1) ⊻ (θo > θs) : false ⊻ (θs > π/2)
        νθfalse =  isincone ? (n%2==1) ⊻ (θo > θs) : true ⊻ (θs > π/2)
        for _i in 1:(2steps)
            α = αvals[_i]
            for _j in 1:(2steps)
                β = βvals[_j]
                rstemp, νrtemp, numroots = rs(α, β, θs, θo, a, true, n)
                evalpol = (rstemp >= rh && rstemp != Inf)
                sintemp, costemp =  evalpol ? calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθtrue, rpeak, p1, p2) : (0., 0.)


                rstemp2, νrtemp2, numroots2 = rs(α, β, π-θs, θo, a, true, n)
                evalpol = (rstemp2 >= rh && rstemp2 != Inf)
                sintemp, costemp =  evalpol ? calcPol(α, β, rstemp2, π-θs, θo, a, B, βfluid, νrtemp2, νθtrue, rpeak, p1, p2) : (0., 0.)


                rsvals[_i][_j] = rstemp
                rsvals3[_i][_j] = rstemp2
                rootvals[_i][_j] = numroots
                sinvals[_i][_j] = sintemp * αmax / 20
                cosvals[_i][_j] = costemp * βmax / 20

                rstemp, νrtemp, numroots = rs(α, β, θs, θo, a, false, n)
                evalpol = (rstemp >= rh && rstemp != Inf)
                sintemp, costemp =  evalpol ? calcPol(α, β, rstemp, θs, θo, a, B, βfluid, νrtemp, νθfalse, rpeak, p1, p2) : (0., 0.)

                rstemp2, νrtemp2, numroots2 = rs(α, β, π-θs, θo, a, false, n)
                evalpol = (rstemp2 >= rh && rstemp2 != Inf)
                sintemp, costemp =  evalpol ? calcPol(α, β, rstemp2, π-θs, θo, a, B, βfluid, νrtemp2, νθfalse, rpeak, p1, p2) : (0., 0.)


                rsvals2[_i][_j] += rstemp
                rsvals4[_i][_j] += rstemp2
                rootvals[_i][_j] = numroots
                sinvals[_i][_j] += sintemp * αmax / 20
                cosvals[_i][_j] += costemp * βmax / 20

                geotypevals[_i][_j] = (Gθ(α, β, a, θs, θo, false,n)[2] ? 1 : 0)

            end
        end
        tickx = reduce(vcat,transpose.(sinvals)) 
        ticky = reduce(vcat,transpose.(cosvals))

        arrows!(ax2, .-αvalstemp, βvalstemp, .-tickx, ticky, arrowhead = ' ')
        if abs(cos(θo)) > abs(cos(θs))
            rsvals .+= rsvals2
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/5:rmax)
        else 
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals)),  levels=rmin:(rmax - rmin)/5:rmax)
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals2)), levels=rmin:(rmax - rmin)/5:rmax)
        end

        if abs(cos(θo)) > abs(cos(θs))
            rsvals3 .+= rsvals4
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals3)),  levels=rmin:(rmax - rmin)/5:rmax)
        else 
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals3)),  levels=rmin:(rmax - rmin)/5:rmax)
            contour!(ax1, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rsvals4)), levels=rmin:(rmax - rmin)/5:rmax)
        end

        heatmap!(ax0, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(rootvals)))
        heatmap!(ax3, .-αvalstemp, βvalstemp, reduce(vcat,transpose.(geotypevals)))

    end
    display(f)
end

a = 1
αmax = 10
βmax = 10
steps = 30
n = 0
θo = (180-60)*π/180
θs = (90)*π/180
back = true

βv = 0.71
θz = π/2
ϕz = -90*π/180
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
    plot(αmax, βmax, rmin, rmax, n, a,steps, θs, θo, B, βfluid, count)
    count +=1
end
a = -0.9986593485349988
θo = 0.229669959563374
θs = 9666474164045953
pa = 3.1440968196839307
βv = 0.9860658551413877
χ = -1.956050891045084
ι = 1.390577214486354
rpeak = 3.5#2.9592973448523763
p1 = 5.5593933604673514
p2 = 2.6630142060526514
n=1
plot(αmax, βmax, rmin, rmax, n, a,steps, θs, θo, B, βfluid, rpeak, p1, p2, pa)



