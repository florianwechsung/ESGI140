# import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib.animation as animation

# def main():
#     numframes = 100
#     numpoints = 10
#     color_data = np.random.random((numframes, numpoints))
#     x, y, c = np.random.random((3, numpoints))

#     fig = plt.figure()
#     scat = plt.scatter(x, y, c=c, s=100)

#     ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes),
#                                   fargs=(color_data, scat))
#     ani.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#     plt.show()

# def update_plot(i, data, scat):
#     scat.set_array(data[i])
#     return scat,

# main()

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def make_runner_video(res_t, res_x, res_rho, width):

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    xmin = np.min(res_x)
    xmax = np.max(res_x)
    xroad = list(range(int(xmin), int(xmax)))
    w0 = width(0)
    plt.plot(xroad, [width(xx) + (w0-width(xx))/2 for xx in xroad])
    plt.plot(xroad, [(w0-width(xx))/2 for xx in xroad])
    scat = plt.scatter([], [], lw=2)
    plt.xlim((xmin, xmax))
    plt.ylim((-1, w0+1))
    yref = np.random.rand(res_x.shape[1])
    line, = plt.plot([], [])

    # initialization function: plot the background of each frame
    def init():
        scat.set_offsets(np.zeros((0, 2)))
        line.set_data([], [])
        return scat,

    # animation function.  This is called sequentially
    def animate(i):
        i = i*3
        print(i)
        x = res_x[i, :]
        y = np.zeros(x.shape)
        for j in range(len(x)):
            w = width(x[j])
            y[j] = yref[j] * w + (w0-w)/2
        z = np.column_stack((x, y))
        # from IPython import embed
        # embed()
        scat.set_offsets(z)
        line.set_data(x, res_rho[i, :])
        return scat,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(res_t)//3, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
