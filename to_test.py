import matplotlib, sys
import matplotlib.pyplot as plt

if __name__ == "__main__" :
    x = sys.argv[1]
    f = open("plot_answer.txt", "r")
    f1 = open("generate_test2.txt", "r")
    
    p=int(f.readline())
    plt.figure()
    for ip in range(p):
        n = int(f.readline())
        n = n + 1
        coord = []
        for i in range(n) :
            co = f.readline().split(' ')
            coord.append([float(co[0]), float(co[1])])

        coord.append(coord[0]) #repeat the first point to create a 'closed loop'
        # coord.reverse()
        xs, ys = zip(*coord) #create lists of x and y values

        # if ip == 0:
        #     plt.plot(xs,ys,linewidth=2, color = 'black') 
        # else:
        plt.fill(xs,ys, linewidth=0.75, edgecolor = 'black') 
    
    
    
    p1 = int(f1.readline())
    coord = []
    for ip2 in range(p1):
        co = f1.readline().split(' ')
        coord.append([float(co[0]), float(co[1])])
    coord.append(coord[0])
    xs2, ys2 = zip(*coord)
    
    plt.plot(xs2, ys2, linewidth = 2, color = 'black')

    font1 = {'color':'purple','size':20}
    plt.title("Value of N = " + x, fontdict = font1, style = 'italic')
    plt.savefig("Outputfiles/output"+sys.argv[2]+".png")
    # plt.show() 
