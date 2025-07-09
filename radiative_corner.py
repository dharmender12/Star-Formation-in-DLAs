# def make_corner_plots( my_chains_matrix ):

#         import numpy as np
#         import pylab as plt

#         N_dim = 6

#         ax_list = []

#         label_list = [ 'log V' , 'log N' , 'log ta' , 'z' , 'log EW', 'Wi'  ]
#         # label_list = [ 'log V' , 'log N' , 'log ta' , 'z' , 'Wi'  ]


#         MAIN_VALUE_mean   = np.zeros(N_dim)
#         MAIN_VALUE_median = np.zeros(N_dim)
#         MAIN_VALUE_MAX    = np.zeros(N_dim)

#         for i in range( 0 , N_dim ):

#             x_prop = my_chains_matrix[ : , i ]

#             x_prop_min = np.percentile( x_prop , 10 )
#             x_prop_50  = np.percentile( x_prop , 50 )
#             x_prop_max = np.percentile( x_prop , 90 )

#             x_min = x_prop_50 - ( x_prop_max - x_prop_min ) * 1.00
#             x_max = x_prop_50 + ( x_prop_max - x_prop_min ) * 1.00

#             mamamask = ( x_prop > x_min ) * ( x_prop < x_max )

#             MAIN_VALUE_mean[  i] = np.mean(       x_prop[ mamamask ] )
#             MAIN_VALUE_median[i] = np.percentile( x_prop[ mamamask ] , 50 )

#             HH , edges_HH = np.histogram( x_prop[ mamamask ] , 30 , range=[ x_prop_min , x_prop_max ] )

#         plt.figure( figsize=(15,15) )

#         Q_top = 80
#         Q_low = 20

#         for i in range( 0 , N_dim ):

#             y_prop = my_chains_matrix[ : , i ]

#             y_prop_min = np.percentile( y_prop , Q_low )
#             y_prop_50  = np.percentile( y_prop , 50 )
#             y_prop_max = np.percentile( y_prop , Q_top  )

#             mask_y = ( y_prop > y_prop_min ) * ( y_prop < y_prop_max )

#             y_min = y_prop_50 - np.std( y_prop[ mask_y ] )
#             y_max = y_prop_50 + np.std( y_prop[ mask_y ] )

#             for j in range( 0 , N_dim ):

#                 if i < j : continue

#                 x_prop = my_chains_matrix[ : , j ]

#                 x_prop_min = np.percentile( x_prop , Q_low )
#                 x_prop_50  = np.percentile( x_prop , 50 )
#                 x_prop_max = np.percentile( x_prop , Q_top )

#                 mask_x = ( x_prop > x_prop_min ) * ( x_prop < x_prop_max )

#                 x_min = x_prop_50 - np.std( x_prop[ mask_x ] )
#                 x_max = x_prop_50 + np.std( x_prop[ mask_x ] )

#                 ax = plt.subplot2grid( ( N_dim , N_dim ) , (i, j)  )

#                 ax_list += [ ax ]

#                 DDX = x_max - x_min
#                 DDY = y_max - y_min

#                 if i==j :

#                     H , edges = np.histogram( x_prop , 30 , range=[x_min,x_max] )

#                     ax.hist( x_prop , 30 , range=[x_min,x_max] , color='cornflowerblue' )

#                     ax.plot( [ MAIN_VALUE_median[i] , MAIN_VALUE_median[i] ] , [ 0.0 , 1e10 ] , 'k--' , lw=2 )

#                     ax.set_ylim( 0 , 1.1 * np.amax(H) )

#                 else :

#                     XX_min = x_min - DDX * 0.2
#                     XX_max = x_max + DDX * 0.2

#                     YY_min = y_min - DDY * 0.2
#                     YY_max = y_max + DDY * 0.2

#                     H , edges_y , edges_x = np.histogram2d( x_prop , y_prop , 30 , range=[[XX_min , XX_max],[YY_min , YY_max]] )

#                     y_centers = 0.5 * ( edges_y[1:] + edges_y[:-1] )
#                     x_centers = 0.5 * ( edges_x[1:] + edges_x[:-1] )

#                     H_min = np.amin( H )
#                     H_max = np.amax( H )

#                     N_bins = 10000

#                     H_Arr = np.linspace( H_min , H_max , N_bins )[::-1]

#                     fact_up_Arr = np.zeros( N_bins )

#                     TOTAL_H = np.sum( H )

#                     for iii in range( 0 , N_bins ):

#                         mask = H > H_Arr[iii]

#                         fact_up_Arr[iii] = np.sum( H[ mask ] ) / TOTAL_H

#                     H_value_68 = np.interp( 0.680 , fact_up_Arr , H_Arr )
#                     H_value_95 = np.interp( 0.950 , fact_up_Arr , H_Arr )

#                     ax.pcolormesh( edges_y , edges_x , H.T , cmap='Blues' )

#                     ax.contour( y_centers, x_centers , H.T , colors='k' , levels=[ H_value_95 ] )
#                     ax.contour( y_centers, x_centers , H.T , colors='r' , levels=[ H_value_68 ] )

#                     X_VALUE =  MAIN_VALUE_median[j]
#                     Y_VALUE =  MAIN_VALUE_median[i]

#                     ax.plot( [ X_VALUE , X_VALUE ] , [    -100 ,     100 ] , 'k--' , lw=2 )
#                     ax.plot( [    -100 ,     100 ] , [ Y_VALUE , Y_VALUE ] , 'k--' , lw=2 )

#                     ax.set_ylim( y_min-0.05*DDY , y_max+0.05*DDY )

#                 ax.set_xlim( x_min-0.05*DDX , x_max+0.05*DDX )

#                 if i==N_dim-1:
#                     ax.set_xlabel( label_list[j] , size=20 )

#                 if j==0 and i!=0 :
#                     ax.set_ylabel( label_list[i] , size=20 )

#                 if j!=0:
#                     plt.setp( ax.get_yticklabels(), visible=False)

#                 if j==0 and i==0:
#                     plt.setp( ax.get_yticklabels(), visible=False)

#                 if i!=len( label_list)-1 :
#                     plt.setp( ax.get_xticklabels(), visible=False)

#         plt.subplots_adjust( left = 0.09 , bottom = 0.15 , right = 0.98 , top = 0.99 , wspace=0., hspace=0.)

#         return None



def make_corner_plots(my_chains_matrix):
    import numpy as np
    import pylab as plt

    N_dim = 5  # Adjusted to 5 since we're removing log_EW

    ax_list = []

    # Updated label list without 'log EW'
    label_list = ['log V', 'log N', 'log $\\tau_\\alpha$', 'z', 'Wi']

    MAIN_VALUE_mean = np.zeros(N_dim)
    MAIN_VALUE_median = np.zeros(N_dim)
    MAIN_VALUE_MAX = np.zeros(N_dim)

    # Loop through the dimensions, skipping the index for log_EW
    for i in range(N_dim):
        x_prop = my_chains_matrix[:, i]

        # x_prop_min = np.percentile(x_prop, 10)
        # x_prop_50 = np.percentile(x_prop, 50)
        # x_prop_max = np.percentile(x_prop, 90)
        x_prop_min = np.percentile(x_prop, 2)
        x_prop_50 = np.percentile(x_prop, 50)
        x_prop_max = np.percentile(x_prop, 98)

        x_min = x_prop_50 - (x_prop_max - x_prop_min) * 1.00
        x_max = x_prop_50 + (x_prop_max - x_prop_min) * 1.00

        mamamask = (x_prop > x_min) & (x_prop < x_max)

        MAIN_VALUE_mean[i] = np.mean(x_prop[mamamask])
        MAIN_VALUE_median[i] = np.percentile(x_prop[mamamask], 50)

        HH, edges_HH = np.histogram(x_prop[mamamask], 30, range=[x_prop_min, x_prop_max])

    plt.figure(figsize=(15, 15))

    Q_top = 80
    Q_low = 20

    for i in range(N_dim):
        y_prop = my_chains_matrix[:, i]

        y_prop_min = np.percentile(y_prop, Q_low)
        y_prop_50 = np.percentile(y_prop, 50)
        y_prop_max = np.percentile(y_prop, Q_top)

        mask_y = (y_prop > y_prop_min) & (y_prop < y_prop_max)

        y_min = y_prop_50 - np.std(y_prop[mask_y])
        y_max = y_prop_50 + np.std(y_prop[mask_y])

        for j in range(N_dim):
            if i < j:
                continue

            x_prop = my_chains_matrix[:, j]

            x_prop_min = np.percentile(x_prop, Q_low)
            x_prop_50 = np.percentile(x_prop, 50)
            x_prop_max = np.percentile(x_prop, Q_top)

            mask_x = (x_prop > x_prop_min) & (x_prop < x_prop_max)

            x_min = x_prop_50 - np.std(x_prop[mask_x])
            x_max = x_prop_50 + np.std(x_prop[mask_x])

            ax = plt.subplot2grid((N_dim, N_dim), (i, j))

            ax_list += [ax]

            DDX = x_max - x_min
            DDY = y_max - y_min

            if i == j:
                H, edges = np.histogram(x_prop, 30, range=[x_min, x_max])

                ax.hist(x_prop, 30, range=[x_min, x_max], color='cornflowerblue')

                ax.plot([MAIN_VALUE_median[i], MAIN_VALUE_median[i]], [0.0, 1e10], 'k--', lw=2)

                ax.set_ylim(0, 1.1 * np.amax(H))

            else:
                XX_min = x_min - DDX * 0.2
                XX_max = x_max + DDX * 0.2

                YY_min = y_min - DDY * 0.2
                YY_max = y_max + DDY * 0.2

                H, edges_y, edges_x = np.histogram2d(x_prop, y_prop, 30, range=[[XX_min, XX_max], [YY_min, YY_max]])

                y_centers = 0.5 * (edges_y[1:] + edges_y[:-1])
                x_centers = 0.5 * (edges_x[1:] + edges_x[:-1])

                H_min = np.amin(H)
                H_max = np.amax(H)

                N_bins = 10000

                H_Arr = np.linspace(H_min, H_max, N_bins)[::-1]

                fact_up_Arr = np.zeros(N_bins)

                TOTAL_H = np.sum(H)

                for iii in range(N_bins):
                    mask = H > H_Arr[iii]
                    fact_up_Arr[iii] = np.sum(H[mask]) / TOTAL_H

                H_value_68 = np.interp(0.680, fact_up_Arr, H_Arr)
                H_value_95 = np.interp(0.950, fact_up_Arr, H_Arr)

                ax.pcolormesh(edges_y, edges_x, H.T, cmap='Blues')

                ax.contour(y_centers, x_centers, H.T, colors='k', levels=[H_value_95])
                ax.contour(y_centers, x_centers, H.T, colors='r', levels=[H_value_68])

                X_VALUE = MAIN_VALUE_median[j]
                Y_VALUE = MAIN_VALUE_median[i]

                ax.plot([X_VALUE, X_VALUE], [-100, 100], 'k--', lw=2)
                ax.plot([-100, 100], [Y_VALUE, Y_VALUE], 'k--', lw=2)

                ax.set_ylim(y_min - 0.05 * DDY, y_max + 0.05 * DDY)

            ax.set_xlim(x_min - 0.05 * DDX, x_max + 0.05 * DDX)

            if i == N_dim - 1:
                ax.set_xlabel(label_list[j], size=20)

            if j == 0 and i != 0:
                ax.set_ylabel(label_list[i], size=20)

            if j != 0:
                plt.setp(ax.get_yticklabels(), visible=False)

            if j == 0 and i == 0:
                plt.setp(ax.get_yticklabels(), visible=False)

            if i != len(label_list) - 1:
                plt.setp(ax.get_xticklabels(), visible=False)

    plt.subplots_adjust(left=0.09, bottom=0.15, right=0.98, top=0.99, wspace=0., hspace=0.)

    return None
