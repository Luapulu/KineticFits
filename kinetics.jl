### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ c769e6fe-59c2-11ec-1d28-23d959d24490
using CairoMakie

# ╔═╡ 91f6a6f4-274b-4f6f-81bc-e0313f78cdca
using XLSX

# ╔═╡ 1e4d24b2-a068-4d29-9061-c7196a2d1d1d
using OrdinaryDiffEq

# ╔═╡ a4ae4581-e7c7-4abb-a844-0e0bb34de766
using LsqFit

# ╔═╡ 4d2a70e0-6852-4125-a271-aba216ac5ed8
using DiffEqParamEstim

# ╔═╡ 2e61a3a8-8bb8-488f-be17-52e361a6c198
using BlackBoxOptim

# ╔═╡ 2d14ff52-5215-499f-aa4f-b4bcfaf1581f
const R = 8.31446261815324

# ╔═╡ a08e9b86-2d75-4eb0-b21e-f8e7fa0e45b8
getdatapaths(dirpath) = filter(readdir(dirpath; join=true)) do p
	return basename(p)[1:4] == "MK01" && split(p, '.')[end] == "xlsx"
end

# ╔═╡ 06c9e5d6-2e60-4713-8362-29761a55b929
datapaths = getdatapaths("../data")

# ╔═╡ 885005f0-60c7-4c72-b14c-2eafb4afb90a
function shiftzeros!(data)
	T = eltype(data)
	for i in eachindex(data)
		data[i] = max(eps(T), data[i])
	end

	return data
end

# ╔═╡ bc6c7f12-484a-4001-9b05-1a6f6c27efd5
function getdata(filepath)
    rawdata = XLSX.readtable(filepath, "Sheet1")

	tidx = findfirst(c -> c == :Time, rawdata[2])
	t = rawdata[1][tidx] |> Vector{Float64}
	t .-= first(t)  # Initial time set to 0

	Tidx = findfirst(c -> c == :Temperature_Cal, rawdata[2])
	T = rawdata[1][Tidx] |> Vector{Float64}
	T .+= 273.15  # to Kelvin

	αidx = findfirst(c -> c == :Conversion_DSC, rawdata[2])
	α = rawdata[1][αidx] |> Vector{Float64}
    return (T = T, α = α, t = t)
end

# ╔═╡ 05d49511-9430-494e-9713-336acef62525
filepath = datapaths[1]

# ╔═╡ 831a4e12-e2a8-4887-9d1c-1fa933adec92
pr=0.5

# ╔═╡ 8ba7c5cc-47ad-40b1-91e8-8bb4562b2f90
data = getdata(filepath)

# ╔═╡ 87d7f36d-c0f5-4e6e-99a6-822e72747dfd
function stripdata(data)
	idx1 = findfirst(x -> x > zero(x), data.α)
	idx2 = findlast(x -> x < one(x), data.α)

	(T = data.T[idx1:idx2], α = data.α[idx1:idx2], t = data.t[idx1:idx2])
end

# ╔═╡ e5d0e591-4bc2-4424-b30f-fe812a6f8a49
let fig = Figure()
	ax = Axis(fig[1, 1])

	lines!(ax, data.T, data.α)

	fig
end

# ╔═╡ 81f11d6e-b418-4b04-8653-4baea77a2087
function linfit(xs, ys)
    linmodel(x, p) = x .* p[1] .+ p[2]
	fit = curve_fit(linmodel, xs, ys, [1.0, 1.0])
	return (a = fit.param[1], b = fit.param[2], fit = fit)
end

# ╔═╡ b4c4af91-de7e-4a6a-bf20-54531ecd46f3
β, T0, Tfit = linfit(data.t, data.T)

# ╔═╡ d9b2b78c-052b-429f-a3b1-061eeeb7fd85
β * 60, T0 - 273.15

# ╔═╡ 216e2df2-8159-46b2-aca1-1b0eeb2b2eba
let fig = Figure()
	ax = Axis(fig[1, 1])

	lines!(ax, data.t, data.T, label="experiment")
	lines!(ax, data.t, β .* data.t .+ T0; label="fit")
	# lines!(ax, data.t, 4 / 60 .* data.t .+ (817.7268 + 273.15))

	axislegend(ax)
	fig
end

# ╔═╡ ac10a1f0-a6d2-4a35-9c5e-0e0e0e12634f
abstract type AbstractModel end

# ╔═╡ 417900e2-9fa0-483d-826a-e74d0b419b41
struct FxModel <: AbstractModel end

# ╔═╡ 58895114-1b00-41cf-a65c-6390e1c3d11d
function makeratefunction(β, T0, pr, ::FxModel)
	return function rate(α, p, t)
		A, B, Ea, ΔH, ΔS, C = p
		T = β * t + T0

		r1 = exp(A - Ea / R / T)
		r2 = 1 - pr * exp((ΔH / T - ΔS) / R)
		r3 = 1 - α
	return r1 * max(r2, 0)^B * max(r3, 0)^C
	end
end

# ╔═╡ 6c3e9753-a624-4d3f-a303-e21a45a39f9d
struct F1Model <: AbstractModel end

# ╔═╡ 035c8b21-0dbb-4652-b598-87ea8711d268
αrate(α, ::F1Model) = max(1 - α, 0)

# ╔═╡ b1cadb47-64d4-4d48-b9dc-e44d2cd7c43b
struct AE3Model <: AbstractModel end

# ╔═╡ f4a948fc-02ac-4111-8ed8-b0342d58a83f
αrate(α, ::AE3Model) = max(3 * (1 - α) * (-log(clamp(1 - α, eps(), 1 - eps())))^(2//3), 0)

# ╔═╡ 6eccb017-8abe-4c88-8afe-3d249f0ff0c0
struct PTModel <: AbstractModel end

# ╔═╡ 89723b6c-8ec7-4a16-9d1a-b93c1303f38e
αrate(α, ::PTModel) = (1 - α) * α

# ╔═╡ eb163472-c8b0-48af-95a0-cb9f6552fb48
struct DJ3Model <: AbstractModel end

# ╔═╡ 8ccb1799-c7f0-497d-8f67-bb3dc285839a
αrate(α, ::DJ3Model) = (c = cbrt(1 - α); 3/2 * c^2 / (1 - c))

# ╔═╡ 98d00c87-099a-4a3c-a5cf-e1300261f29b
struct GBModel <: AbstractModel end

# ╔═╡ 55760264-1ca1-43ed-a0ee-506d39098962
αrate(α, ::GBModel) = (c = cbrt(1 - α); 3/2 * c / (1 - c))

# ╔═╡ d73d9a7f-1f05-4355-9652-6b33bb27c78a
function makeratefunction(β, T0, pr, model::AbstractModel)
	return function rate(α, p, t)
		A, B, Ea, ΔH, ΔS = p
		T = β * t + T0

		r1 = exp(A - Ea / R / T)
		r2 = 1 - pr * exp((ΔH / T - ΔS) / R)
		r3 = αrate(α, model)
	return r1 * max(r2, 0)^B * r3
	end
end

# ╔═╡ 536a0f3d-70c3-4e7d-a136-d1262f491916
const MODELS = (FxModel(), F1Model(), AE3Model(), PTModel(), DJ3Model(), GBModel())

# ╔═╡ dbc79c4e-5ca2-47df-a7e4-19dc35149a18
makeratefunction(4, 500, 0.3, FxModel())(0.5, [20, 5, 2e5, 1.7e5, 150, 1.2], 456)

# ╔═╡ edab6163-e0f3-40bc-af8e-1d071ce9ec54
function fitkinetics(rate, ts, αs, p0, bounds, alg; losskwargs=(), solkwargs=(), optkwargs=())
	prob = ODEProblem(rate, first(αs), (first(ts), last(ts)), p0)
	loss = L2Loss(ts, αs)
	loss_objective = build_loss_objective(prob, alg, loss; solkwargs...)
	bboptimize(loss_objective, p0; SearchRange=bounds, optkwargs...)
end

# ╔═╡ 981584bd-3a65-4f99-83be-0095ffd61ea6
rate = makeratefunction(β, T0, pr, AE3Model())

# ╔═╡ fc1f8150-0b88-4d29-abb8-fceeb008bca0
p0 = [67.75974468855404, 10.891425746932146, 696081.8832423445, 131682.74431826902, 264.8105989317905]

# ╔═╡ 588951eb-d96e-4587-8eba-d78527aa2f97
let fig = Figure()
	ax = Axis(fig[1, 1])
	# ylims!(ax, -0.01, 0.015)

	dαs = [0; (data.α[2:end] .- data.α[1:end-1]) ./
		(data.t[2:end] .- data.t[1:end-1])]
	r0s = map(zip(data.α, data.t)) do (a, t)
		rate(a, p0, t)
	end

	lines!(ax, data.t, dαs; label="experiment")
	lines!(ax, data.t, r0s; label="model 0")

	axislegend(ax)

	fig
end

# ╔═╡ f432693b-0c8b-4f23-b4dc-99e2707fb9f8
res = fitkinetics(rate, data.t, data.α,
	p0, [(-50, 200), (-50, 50), (1e4, 2e6), (1e5, 2e5), (50, 300)],
	BS3(); losskwargs=(differ_weight=0,), solkwargs=(maxiters=1000, reltol=2e-3),
	optkwargs=(MaxTime=1, Method=:adaptive_de_rand_1_bin_radiuslimited)
)

# ╔═╡ 43f37494-6247-4b65-bdd4-2f5983981d9d
# C bounds = (0.0, 5.0)

# ╔═╡ b069488b-e485-46e1-966d-c1308514c4fa
bc = best_candidate(res)

# ╔═╡ e912818f-ee76-490d-934d-cd8cc2f52878
string(bc)

# ╔═╡ c6740e72-3870-4444-bdac-9d7529d37b36
best_fitness(res)

# ╔═╡ 2d4b3cbb-1f81-4962-bc64-55bb75e929aa
tspan = (first(data.t), last(data.t))

# ╔═╡ 58ad0dbe-94c7-4d31-b9fe-2c9ef99377d9
prob0 = ODEProblem(rate, first(data.α), tspan, p0)

# ╔═╡ 93fdc768-03cd-45ae-a1a7-69cfba984fdc
sol0 = solve(prob0, Vern7(); saveat=data.t,save_everystep=false,dense=false, reltol=1e-4);

# ╔═╡ db6e32dd-6c8a-4d01-a869-98d3d8d9a5a6
sum(abs2, sol0.(data.t) .- data.α)

# ╔═╡ 2b958eb1-e7d7-4abf-8537-34f5769d5482
bestprob = ODEProblem(rate, first(data.α), tspan, bc)

# ╔═╡ 1fa0c8a6-9e41-4840-a17b-8727445156b1
bestsol = solve(bestprob, Vern7(); saveat=data.t,save_everystep=false,dense=false, reltol=1e-4);

# ╔═╡ 7d5219c6-ad93-4dea-87e3-a092f2288e63
sum(abs2, bestsol.(data.t) .- data.α)

# ╔═╡ 7c0fe222-7413-4187-9bf2-e425f3af473c
let fig = Figure()
	ax = Axis(fig[1, 1])
	xlims!(ax, tspan)
	ylims!(ax, (-0.05, 1.05))

	xs = data.t
	ys = data.α

	lines!(ax, xs, ys; label="experiment")
	lines!(ax, xs, sol0.(xs); label="model p0")
	lines!(ax, xs, bestsol.(xs); label="best model")

	axislegend(ax)

	fig
end

# ╔═╡ 7fb59c6b-80cc-4e52-9cad-6a6ac8a635b9
let fig = Figure()
	ax = Axis(fig[1, 1])
	# ylims!(ax, -0.01, 0.015)

	dαs = [0; (data.α[2:end] .- data.α[1:end-1]) ./
		(data.t[2:end] .- data.t[1:end-1])]
	r0s = map(zip(data.α, data.t)) do (a, t)
		rate(a, p0, t)
	end
	brs = map(zip(data.α, data.t)) do (a, t)
		rate(a, bc, t)
	end

	lines!(ax, data.t, dαs; label="experiment")
	lines!(ax, data.t, r0s; label="model 0")
	lines!(ax, data.t, brs; label="best model")

	axislegend(ax)

	fig
end

# ╔═╡ 8dca68e6-b22a-47f0-b5c2-59936dbf40bd
let fig = Figure()
	ax = Axis(fig[1, 1])
	xlims!(ax, tspan)
	ylims!(ax, (-0.1, 0.1))

	xs = data.t
	ys = data.α

	lines!(ax, xs, zeros(length(xs)))
	lines!(ax, xs, sol0.(xs) .- ys)
	lines!(ax, xs, bestsol.(xs) .- ys)

	fig
end

# ╔═╡ 11045900-4a1f-430b-8412-f8dd63508818
const RATES = makeratefunction.(β, T0, 0.5, MODELS)

# ╔═╡ d7bf524a-0092-4775-bd9d-623e3a59019a
modeltostring(model::AbstractModel) = split(string(model), '.')[end]

# ╔═╡ 78d7ce10-9218-4a00-acce-40755a97328d
let fig = Figure()
	ax = Axis(fig[1, 1])
	ylims!(ax, -0.01, 0.02)

	lines!(ax, data.t, dαs; label="experiment")

	for model in MODELS
		rf = makeratefunction(β, T0, 0.5, model)
		rs = map(zip(data.α, data.t)) do (a, t)
			rf(a, p0, t)
		end
		lines!(ax, data.t, rs; label=modeltostring(model))
	end

	axislegend(ax)

	fig
end

# ╔═╡ Cell order:
# ╠═c769e6fe-59c2-11ec-1d28-23d959d24490
# ╠═91f6a6f4-274b-4f6f-81bc-e0313f78cdca
# ╠═1e4d24b2-a068-4d29-9061-c7196a2d1d1d
# ╠═a4ae4581-e7c7-4abb-a844-0e0bb34de766
# ╠═4d2a70e0-6852-4125-a271-aba216ac5ed8
# ╠═2e61a3a8-8bb8-488f-be17-52e361a6c198
# ╠═2d14ff52-5215-499f-aa4f-b4bcfaf1581f
# ╠═a08e9b86-2d75-4eb0-b21e-f8e7fa0e45b8
# ╠═06c9e5d6-2e60-4713-8362-29761a55b929
# ╠═885005f0-60c7-4c72-b14c-2eafb4afb90a
# ╠═bc6c7f12-484a-4001-9b05-1a6f6c27efd5
# ╠═05d49511-9430-494e-9713-336acef62525
# ╠═831a4e12-e2a8-4887-9d1c-1fa933adec92
# ╠═8ba7c5cc-47ad-40b1-91e8-8bb4562b2f90
# ╠═87d7f36d-c0f5-4e6e-99a6-822e72747dfd
# ╠═e5d0e591-4bc2-4424-b30f-fe812a6f8a49
# ╠═81f11d6e-b418-4b04-8653-4baea77a2087
# ╠═b4c4af91-de7e-4a6a-bf20-54531ecd46f3
# ╠═d9b2b78c-052b-429f-a3b1-061eeeb7fd85
# ╠═216e2df2-8159-46b2-aca1-1b0eeb2b2eba
# ╠═ac10a1f0-a6d2-4a35-9c5e-0e0e0e12634f
# ╠═417900e2-9fa0-483d-826a-e74d0b419b41
# ╠═58895114-1b00-41cf-a65c-6390e1c3d11d
# ╠═d73d9a7f-1f05-4355-9652-6b33bb27c78a
# ╠═6c3e9753-a624-4d3f-a303-e21a45a39f9d
# ╠═035c8b21-0dbb-4652-b598-87ea8711d268
# ╠═b1cadb47-64d4-4d48-b9dc-e44d2cd7c43b
# ╠═f4a948fc-02ac-4111-8ed8-b0342d58a83f
# ╠═6eccb017-8abe-4c88-8afe-3d249f0ff0c0
# ╠═89723b6c-8ec7-4a16-9d1a-b93c1303f38e
# ╠═eb163472-c8b0-48af-95a0-cb9f6552fb48
# ╠═8ccb1799-c7f0-497d-8f67-bb3dc285839a
# ╠═98d00c87-099a-4a3c-a5cf-e1300261f29b
# ╠═55760264-1ca1-43ed-a0ee-506d39098962
# ╠═536a0f3d-70c3-4e7d-a136-d1262f491916
# ╠═dbc79c4e-5ca2-47df-a7e4-19dc35149a18
# ╠═edab6163-e0f3-40bc-af8e-1d071ce9ec54
# ╠═981584bd-3a65-4f99-83be-0095ffd61ea6
# ╠═fc1f8150-0b88-4d29-abb8-fceeb008bca0
# ╟─588951eb-d96e-4587-8eba-d78527aa2f97
# ╠═f432693b-0c8b-4f23-b4dc-99e2707fb9f8
# ╠═43f37494-6247-4b65-bdd4-2f5983981d9d
# ╠═b069488b-e485-46e1-966d-c1308514c4fa
# ╠═e912818f-ee76-490d-934d-cd8cc2f52878
# ╠═c6740e72-3870-4444-bdac-9d7529d37b36
# ╠═2d4b3cbb-1f81-4962-bc64-55bb75e929aa
# ╠═58ad0dbe-94c7-4d31-b9fe-2c9ef99377d9
# ╠═93fdc768-03cd-45ae-a1a7-69cfba984fdc
# ╠═db6e32dd-6c8a-4d01-a869-98d3d8d9a5a6
# ╠═2b958eb1-e7d7-4abf-8537-34f5769d5482
# ╠═1fa0c8a6-9e41-4840-a17b-8727445156b1
# ╠═7d5219c6-ad93-4dea-87e3-a092f2288e63
# ╠═7c0fe222-7413-4187-9bf2-e425f3af473c
# ╟─7fb59c6b-80cc-4e52-9cad-6a6ac8a635b9
# ╠═8dca68e6-b22a-47f0-b5c2-59936dbf40bd
# ╠═11045900-4a1f-430b-8412-f8dd63508818
# ╠═d7bf524a-0092-4775-bd9d-623e3a59019a
# ╠═78d7ce10-9218-4a00-acce-40755a97328d
