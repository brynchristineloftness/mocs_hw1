import matplotlib.pyplot as plt
import numpy as np

def euler(t,h,x,f,f_args):

    """
    given a function f = ,
    return x(t+h), given by x(t)+h*f(x,t)

    args:
        t: current time
        x: current x(t)
        f: a function x'(t) as f(x,t)
        h: size of step
    """

    return x + h*f(x,t,*f_args)


def heun(t,h,x,f,f_args):
    
    """
    given a function f = ,
    return and approximation for x(t+h),
         given by x(t)+h/2*( f(x,t) + euler(t,h,x,f,f_args) )

    args:
        t: current time
        x: current x(t)
        f: a function x'(t) as f(x,t)
        h: size of step
    """

    euler_x = f(x,t,*f_args)
    return x + (h/2)*(euler_x + f(euler_x,t+h,*f_args))


def delta_I(I,t,beta,N=100,gamma=0.25):

    S = N-I
    return beta*I*S - gamma*I


def generate_data(integration_func,betas,time_steps):

    """
    Given an integration technique and a series of model
    parameters, run the SIS model for 50 time steps.

    Returns [time, I_t] as an array of number of infections at time t
    """

    data = []

    for beta in betas:
        for h in time_steps:
            
            #specify initial condition
            
            I_t = [I_0]

            #run simulation
            time = h*np.arange(50)
            for t in time:
                I_t.append(integration_func(t,h,I_t[-1],delta_I,f_args=[beta]))

            # add time of final state
            time = np.append(time, [time[-1]+h])

            data.append((time,np.array(I_t)))

    return data

def format_data(data):

    t_data,I_data = map(np.array, zip(*data))
    I_data = np.array(I_data).reshape((3,3,51))
    t_data = np.array(t_data).reshape((3,3,51))

    return t_data,I_data

if __name__ == "__main__":

    N=100
    I_0 = 10

    time_steps = [0.01, 0.5, 2.0]
    betas = [0.03, 0.06, 0.1]
    euler_data = generate_data(euler,betas,time_steps)
    heun_data = generate_data(heun,betas,time_steps)

    euler_t, euler_I = format_data(euler_data)
    heun_t, heun_I = format_data(heun_data)

    fig, axes = plt.subplots(3,3,figsize=(9,9))

    for i in range(3):
        for j in range(3):
            ax = axes[i,j]

            time = euler_t[i][j]
            I_t_euler = euler_I[i][j]
            I_t_heun = heun_I[i][j]

            #S_t = 100- I_t

            ax.plot(time,I_t_euler,label='I(t) Euler')
            ax.plot(time,I_t_heun,label='I(t) Heun')
            ax.set_xlabel('Time (t)')
            ax.set_ylabel('Num People')  
            ax.legend()
            ax.set_ylim(0,100)
            ax.set_title(f"Beta = {betas[i]}, Time Step = {time_steps[j]}",
                fontsize = 10)



    fig.tight_layout()
    plt.savefig('heun_euler_comparison.png')
    plt.show()





            