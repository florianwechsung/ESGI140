import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def make_runner_video(res_t, res_x, res_v, res_rho, width, num_frames, xmax=None):
    fak = len(res_t)//num_frames
    num_frames = len(res_t)//fak
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    xmin = np.min(res_x)
    if xmax is None:
        xmax = np.max(res_x)
    xroad = list(range(int(xmin), int(xmax)))
    w0 = width(0)
    plt.plot(xroad, [width(xx) + (w0-width(xx))/2 for xx in xroad])
    plt.plot(xroad, [(w0-width(xx))/2 for xx in xroad])
    x = res_x[0, :]
    y = np.zeros(x.shape)
    c = np.asarray([[c, 0, 0] for c in res_v[0, :]])[:,0]
    scat = plt.scatter(x, y, c=np.random.rand(len(x)), s=1)
    plt.xlim((xmin, xmax))
    plt.ylim((-1, w0+1))
    yref = np.random.rand(res_x.shape[1])
    line, = plt.plot([], [])

    # initialization function: plot the background of each frame
    def init():
        scat.set_offsets(np.zeros((0, 2)))
        line.set_data([], [])
        return scat, line

    # animation function.  This is called sequentially
    def animate(i):
        i = i*fak
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
        c = np.asarray([[cc, 0, 0] for cc in res_v[i, :]])[:,0]/5.
        scat.set_array(c)
        plt.title("Time " + str(res_t[i]))
        # line.set_data(x, res_rho[i, :])
        return scat, line

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=num_frames, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
