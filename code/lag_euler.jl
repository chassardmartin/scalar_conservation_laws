# Julia Version of the lagrangian-eulerian code from [Abreu et al., 22] 

using Plots 
using LinearAlgebra 
using LaTeXStrings


function compute_up_conv(x, v, w, rho_n, dx, j)
    N = length(x) 
    _w = collect([w(k*dx[j]) for j=-N:1]) 
    to_sum = []

    for k =-N:1
        if j+k >= 0
            rho_term = rho_n[j+k] 
        else
            rho_term = rho_n[N-1+j+k]
        end
        push!(to_sum, _w[k+N-1]*rho_term)
    end
    return v(dx[j] * sum(to_sum))
end

function up_conv(x, v, w, rho_n, dx)
    N = length(x) 
    return collect(
        [compute_up_conv(x, v, w, rho_n, dx, j) for j=1:N-2]
    )
end

function compute_down_conv(x, v, w, rho_n, dx, j)
    N = length(x) 
    _w = collect([w(k*dx[j]) for k=1:N]) 
    to_sum = []

    for k =1:N
        if j+k <= N 
            rho_term = rho_n[j+k] 
        else
            rho_term = rho_n[-(N-1)+j+k]
        end
        push!(to_sum, _w[k]*rho_term)
    end
    return v(dx[j] * sum(to_sum))
end

function down_conv(x, v, w, rho_n, dx)
    N = length(x) 
    return collect(
        [compute_down_conv(x, v, w, rho_n, dx, j) for j=1:N-2]
    )
end

function evolve_space_step(x, dx, dt, V_n, rho_n, f)
    N = length(x)
    new_dx = copy(dx) 
    new_dx[1:N-3] = collect(
        [dx[j] + dt * (f(rho_n[j+1]*V_n[j+1]*rho_n[j+1]^(-1))
                        - f(rho_n[j]*V_n[j]*rho_n[j]^(-1)))
                        for j=1:N-3]
    )
    return new_dx 
end

function locally_update_rho(rho_n, new_dx, dt, V_n, f, j)
    new_rho = (rho_n[j-1] + 2*rho_n[j] + rho_n[j+1]) / 4 
    if (rho_n[j-1]) != 0 && (rho_n[j] != 0) && (rho_n[j+1] != 0)
        new_rho += (dt/4) * (
            (new_dx[j-1]^(-1)) * (f(rho_n[j-1])*V_n[j-1]*rho_n[j-1]^(-1) + f(rho_n[j])*V_n[j]*rho_n[j]^(-1)) * (rho_n[j-1] + rho_n[j])
            - (new_dx[j+1]^(-1)) * (f(rho_n[j])*V_n[j]*rho_n[j]^(-1) + f(rho_n[j+1])*V_n[j+1]*rho_n[j+1]^(-1)) * (rho_n[j] + rho_n[j+1])
        )
    else
        c = 0 
    end
    return new_rho 
end

function update_rho(x, rho_n, new_dx, dt, V_n, f)
    N = length(x)
    new_rho = copy(rho_n)
    new_rho[2:N-3] = collect(
        [locally_update_rho(rho_n, new_dx, dt, V_n, f, j)
         for j=2:N-3]
    )
    ### Absorbing BC : extending the solution constantly on the right boundary 
    last_value = new_rho[N-3]
    new_rho[N-1:end] = collect([last_value, last_value])
    return new_rho
end

########## Solver  ####################

function solve(rho_0, v, w, f, N, x_max, x_min, dt, T, convolution)

    dx_element = (x_max - x_min) / N
    dx_zero = repeat([dx_element], N-1)
    x = collect([x_min + j*dx_element for j=1:N])
    sol = [] 
    convolution_method = convolution

    rho_zero = collect(rho_0.(x))
    if convolution_method == "upstream"
        V_zero = up_conv(x, v, w, rho_zero, dx_zero)
    else
        V_zero = down_conv(x, v, w, rho_zero, dx_zero)
    end

    dx = dx_zero 
    rho = rho_zero 
    V = V_zero 
    n_iterations = floor(Int, (T/dt)) 

    for n = 1:n_iterations 
        
        if convolution_method == "upstream"
            new_V = up_conv(x, v, w, rho, dx)
        else
            new_V = down_conv(x, v, w, rho, dx)
        end

        new_dx = evolve_space_step(x, dx, dt, V, rho, f)
        new_rho = update_rho(x, rho, new_dx, dt, new_V, f)

        if mod(n, n_iterations//10) == 0 
            println("CFL:")
            println("max Vn :", maximum(abs.(new_V)))
            println("RHS :", (1/8) * maximum(new_dx) / dt)
            println("")
        end

        push!(sol, new_rho)
        V = new_V 
        dx = new_dx 
        rho = new_rho 
    end
    return sol
end 


############### Main ##################### 


#### Sedimentation Model #########
# equation = "sedimentation"
# eta = 1.0;
# x_min = 0.0;
# x_max = 4.0;  
# N = 4000; 
# h = (x_max - x_min) / N;
# dt = h / 40; 
# _x = collect([x_min + j*h for j=0:N-1])

# σ = 0.2
# rho_0(x) = exp(-(x-2.0)^2/(2*σ^2));

# function k(x)
#     if (x>-2) && (x<2)
#         return 3/8 * (1 - x^2 /4)
#     else
#         return 0.0
#     end
# end

# w(x) =  k(x/eta) / eta; 
# f(x) = x * (1-x); 
# v(x) = (1-x)^4; 
####################################

#### Burgers #################################
equation = "burgers"
eta = 1.0; 
x_min = -10; 
x_max = 10;
N = 1000; 
h = (x_max - x_min) / N;
dt = h / 40;
_x = collect([x_min + j*h for j=0:N-1]);

w(x) = 1/eta; 
f(x) = 1/2 * x^2;
v(x) = 1 - x;
rho_0(x) = 4 * ((1/(7*pi)^0.5) * exp(- x^2/7));
##############################################

### Analytic Rectangle ##################
# equation = "rectangle";
# eta = 1.0;
# x_min = -1;
# x_max = 2;
# N = 2000;
# h = (x_max - x_min) / N;
# dt = h/40;
# _x = collect([x_min + j*h for j=0:N-1]);

# function w(x)
#     if x>=0 && x<1
#         return 1.0
#     else
#         return 1e-10
#     end
# end

# f(x) = x;
# v(x) = x;

# function rho_0(x)
#     if x>=0 && x<1
#         return 1.0
#     else
#         return 1e-10
#     end
# end
###########################################

### Comparison to Blandin #################
# equation = "blandin"
# eta = 0.1;
# x_min = -1;
# x_max = 1;
# N = 3000;
# h = (x_max - x_min) / N;
# dt = h/40;
# _x = collect([x_min + j*h for j=0:N-1]);

# function w(x)
#     if x>=0 && x<eta
#         return 1/eta
#     else
#         return 1e-10
#     end
# end

# f(x) = x;
# v(x) = 1-x;

# function rho_0(x)
#     if x<0
#         return 0.4
#     else 
#         return 0.9
#     end
# end
##############################################


T = 0.1;
@time solutions = solve(rho_0, v, w, f, N, x_max, x_min, dt, T, "downstream");

plot(_x, solutions[1], label="inital cond", title= L"solution $\rho(x,t)$");
n_curves = 10;
for i=1:(n_curves-1)
    element = length(solutions) // n_curves;
    index = floor(Int, i * element)
    time = string(round(index * dt, digits=3))
    _label = "t=" * time * "s"
    plot!(_x, solutions[index], label = _label)
end
final_label = "T=" * string(T) * "s"
plot!(_x, solutions[end], label=final_label); 
xlabel!(L"$x$");

image_name = "../graphs/julia_example_" * equation * "_T=" * string(T) * "_N=" * string(N) * "_dt=" * string(round(dt, digits=3)) * "solving_time=13s.png";
savefig(image_name);