using CGP
using Logging

# inputs = Dict("cart_pole"=>4, "mountain_car"=>2, "pendulum"=>3)
# outputs = Dict("cart_pole"=>2, "mountain_car"=>1, "pendulum"=>1)

function cart_pole(c::Chromosome)
    gravity = 9.81; m_cart = 1.0; m_pole = 0.1; tau = 0.02; force_mag = 10.0
    theta_threshold = 12 * 2 * pi / 360; x_threshold = 2.4; pole_length = 0.5
    m_total = m_cart + m_pole
    limits = [x_threshold, force_mag, theta_threshold, 2*pi]
    reward = 0
    state = 0.1*rand(4)-0.05
    x, x_dot, theta, theta_dot = state

    while ((-x_threshold <= x <= x_threshold) &&
           (-theta_threshold <= theta <= theta_threshold) &&
           reward < 100)
        # get action
        inps = [-1.0; max.(min.(state ./ limits, 1.0), -1.0)[:]]
        action = indmax(process(c, inps))
        force = force_mag * (2 * (action - 1) - 1)
        reward += 1

        # update
        temp = (force + m_pole * pole_length * theta_dot * theta_dot * sin(theta))/m_total
        theta_acc = ((gravity * sin(theta) - cos(theta) * temp) /
                     (pole_length * (4.0/3 - m_pole * cos(theta) * cos(theta) / m_total)))
        x_acc = temp - m_pole * pole_length * theta_acc * cos(theta) / m_total
        x += tau * x_dot
        x_dot += tau * x_acc
        theta += tau * theta_dot
        theta_dot += tau * theta_acc
        state = [x, x_dot, theta, theta_dot]
        # println(reward, " ", action, " ", state, " ", inps)
    end
    reward / 100
end

function mountain_car(c::Chromosome)
    min_x = -1.2; max_x = 0.6; max_v = 0.07
    goal_position = 0.45; power = 0.0015
    steps = 0
    x = -0.6 + 0.2 * rand(); v = 0.0
    best_x = x
    while ((x <= goal_position) && (steps < 100))
        inps = [0.0, 0.0, 2*((x-min_x)/(max_x-min_x))-1.0, v / max_v, 0.0]
        force = process(c, inps)[1]
        v = min(max(v + force * power - 0.0025 * cos(3 * x), -max_v), max_v)
        x = min(max(x + v, min_x), max_x)
        if (x == min_x) && (v < 0)
            v = 0
        end
        steps += 1
        best_x = max(best_x, x)
        # println(steps, " ", inps, " ", force, " ", x, " ", v, " ")
    end
    min(best_x / goal_position, 1.0)
end

function pendulum(c::Chromosome)
    max_speed = 8; max_torque = 2.0; dt = 0.5;
    g = 10.; m = 1.; l = 1.

    th = 2*pi*rand() - pi
    thdot = 2*rand() - 1.0

    cost = 0.0; steps = 0;
    while steps < 100
        inps = [1.0, cos(th), sin(th), thdot / max_speed, 0.0]
        u = max_torque * process(c, inps)[2]

        angle = (((th+pi) % (2*pi)) - pi)
        cost -= angle^2 + (.1*thdot)^2 + 0.001*(u^2)

        thdot = min(max(thdot + (-3*g/(2*l) * sin(th +pi) + 3 ./(m*l^2)*u) * dt,
                        -max_speed), max_speed)
        th = th + thdot * dt
        steps += 1
        # println(steps, " ", inps, " ", u, " ", cost)
    end
    max(min((cost + 4000) / 3850, 1.0), 0.0)
    # min((-100-cost)/-100, 1.0)
end

function combined(c::Chromosome)
    score = cart_pole(c); CGP.reset!(c)
    score += mountain_car(c); CGP.reset!(c)
    score += pendulum(c)
    score
end

function repeat(c::Chromosome, func::Function)
    score = 0
    for i = 1:10
        score += func(c)
    end
    score/10
end

seed = 0
log = "log"
fitness = "pendulum"
if length(ARGS) > 0; seed = parse(Int64, ARGS[1]); end
if length(ARGS) > 1; log = ARGS[2]; end
if length(ARGS) > 2; fitness = ARGS[3]; end

CGP.Config.init("cfg/base.yaml")
CGP.Config.init("cfg/classic.yaml")
# CGP.Config.init("cfg/test.yaml")

Logging.configure(filename=log, level=INFO)
# nin = inputs[fitness]; nout = outputs[fitness]
nin = 5; nout = 2;
fit = x->repeat(x, eval(parse(fitness)))

for ea in CGP.EAs
    dists = CGP.distances[:]
    crosses = CGP.crossovers[:]
    if ea!=cgpneat
        dists = CGP.distances[1:1]
    end
    if ea!=GA
        crosses = CGP.crossovers[1:1]
    end
    for ct in CGP.CTYPES
        for mut in CGP.mutations
            for cross in crosses
                for dist in dists
                    srand(seed)
                    maxfit, best = ea(ct, nin, nout, fit;
                                      f_mutate=mut, f_crossover=cross,
                                      f_distance=dist, seed=seed)
                end
            end
        end
    end
end
