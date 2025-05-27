using DelimitedFiles, Plots

N_levels = 12
N_systems = 3
data = readdlm("energies_q2.txt")
data_res = reshape(data, N_levels+1, N_systems)
x = data_res[1, :]
y = data_res[2:end, :]

p = scatter()
xlabel!(p, "System size")
ylabel!(p, "Energy")
for i=1:N_levels
    scatter!(x, y[i, :], label="Level: $i")
end
display(p)
savefig(p, "plot_q2.png")