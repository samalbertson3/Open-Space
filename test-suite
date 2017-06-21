def test():
    #Displays orbital elements from calculate(). Compares certain elements to known values.
    with open('test_data.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            for i in range(0,len(row)):
                row[i] = float(row[i])
            print '--------------------------------------------------'
            a,e,i,omega,w,nu = calculate(row[0],row[1],row[2],row[3],row[4],row[5],row[6])
            print 'a = ', a, '(', a/a_ref, ')'
            print 'e = ', abs_val(e), '(', abs_val(e)/e_ref, ')'
            print 'i = ', i
            print 'omega = ', omega
            print 'w = ', w
            print 'nu = ', nu
            print '--------------------------------------------------'

def err_test(r0,lng,lat,v0,theta,phi,t):
    #Determines range of possible outputs of calculate(), given estimated error bounds.
    a_set = []
    e_set = []
    i_set = []
    omega_set = []
    w_set = []
    for r1 in [r0,r0+r_err,r0-r_err]:
        for lng1 in [lng,lng+lng_err,lng-lng_err]:
            for lat1 in [lat,lat+lat_err,lat-lat_err]:
                for v1 in [v0,v0+v_err,v0-v_err]:
                    for theta1 in [theta,theta+theta_err,theta-theta_err]:
                        for phi1 in [phi,phi+phi_err,phi-phi_err]:
                            a,e,i,omega,w = calculate(r1,lng1,lat1,v1,theta1,phi1,t)
                            a_set.append(a)
                            e_set.append(abs_val(e))
                            i_set.append(i)
                            omega_set.append(omega)
                            w_set.append(w)
    print r0,rad_to_deg(adj_lng(lng,t)),lat,v0,theta,phi,t
    print '--------------------------------------------------'
    print 'a = ', [int(min(a_set)),int(mean(a_set)),int(max(a_set))]
    print 'e = ', [min(e_set),mean(e_set),max(e_set)]
    print 'i = ', [min(i_set),mean(i_set),max(i_set)]
    print 'omega = ', [min(omega_set),mean(omega_set),max(omega_set)]
    print 'w = ', [min(w_set),mean(w_set),max(w_set)]
    print '--------------------------------------------------'

def test2():
    #Given a set of test inputs, generates a case-by-case estimate of uncertainty in orbital parameters from calculate().
    with open('test_data.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            for i in range(0,len(row)):
                row[i] = float(row[i])
            err_test(row[0],row[1],row[2],row[3],row[4],row[5],row[6])

def test3(r0,lng0,lat0,v0,theta0,phi0,t):
    #Gives values for plot of [arbitrary orbital element] vs. [arbitrary input].
    e_r = []
    r = iterate(r0-500*r_err,1,r0+499*r_err)
    for i in r:
        a,e,i,omega,w = calculate(i,lng0,lat0,v0,theta0,phi0,t)
        e_r.append(omega)
    e_lng = []
    lng = iterate(0,lng_err,2*pi)
    for i in lng:
        a,e,i,omega,w = calculate(r0,i,lat0,v0,theta0,phi0,t)
        e_lng.append(omega)
    e_lat = []
    lat = iterate(0,lat_err,2*pi)
    for i in lat:
        a,e,i,omega,w = calculate(r0,lng0,i,v0,theta0,phi0,t)
        e_lat.append(omega)
    e_v = []
    v = iterate(v0-500*v_err,v_err,v0+499*v_err)
    for i in v:
        a,e,i,omega,w = calculate(r0,lng0,lat0,i,theta0,phi0,t)
        e_v.append(omega)
    e_theta = []
    theta = iterate(0,theta_err,360)
    for i in theta:
        a,e,i,omega,w = calculate(r0,lng0,lat0,v0,i,phi0,t)
        e_theta.append(omega)
    e_phi = []
    phi = iterate(-89,1,89)
    for i in phi:
        a,e,i,omega,w = calculate(r0,lng0,lat0,v0,theta0,i,t)
        e_phi.append(omega)
    with open('test_results.csv','wb') as f:
        w = csv.writer(f)
        w.writerows([phi,e_phi])
        #for i in range(0,len(e_r)):
        #    f.writerows([str(r[i]),str(e_r[i])])
    return phi,e_phi

def test4():
    rs = []
    lngs = []
    lats = []
    vs = []
    thetas = []
    phis = []
    r0 = []
    lng0 = []
    lat0 = []
    v0 = []
    theta0 = []
    phi0 = []
    with open('test_data.csv','rb') as f:
        reader = csv.reader(f)
        for row in reader:
            for i in range(0,len(row)):
                row[i] = float(row[i])
            r0.append(row[0])
            lng0.append(row[1])
            lat0.append(row[2])
            v0.append(row[3])
            theta0.append(row[4])
            phi0.append(row[5])
            a,e,i,omega,w,nu = calculate(row[0],row[1],row[2],row[3],row[4],row[5],row[6])
            r,lng,lat,v,theta,phi = heading(a,e,i,omega,w,nu,row[6])
            rs.append(r)
            lngs.append(lng)
            lats.append(lat)
            vs.append(v)
            thetas.append(theta)
            phis.append(phi)
    return sub(r0,rs),sub(lng0,lngs),sub(lat0,lats),sub(v0,vs),sub(theta0,thetas),sub(phi0,phis)
