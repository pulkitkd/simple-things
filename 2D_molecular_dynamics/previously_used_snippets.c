/*
    for(int i = 0; i < 2000; i++){
    write_position(position, a);
    get_acceleration(position, acceleration, PE);
    velocity = velocity + dt * acceleration;
    position = position + dt * velocity + 0.5 * dt * dt * acceleration;
    boundary_condition(position, L);
    cout<<i<<endl;
    //cout<<"velocity"<<endl<<velocity;
    a++;
    }
    cout<<position;
    cout<<velocity;
    cout<<acceleration;


    while(t < T)
    {
        write_position(position, a);
        get_acceleration(position, acceleration, PE);
        velocity = velocity + dt * acceleration;
        position = position + dt * velocity + 0.5 * dt * dt * acceleration;
        //boundary conditions
        for(int i = 0; i < N; i++)
            for(int j = 0; j < D; j++)
                {
                    if(position(i,j) > 0.99*L)
                    position(i,j) = position(i,j) - L;
                    if(position(i,j) < 0.01)
                    position(i,j) = position(i,j) + L;
                }
//        KE = 0;
//        for(int i = 0; i < N; i++)
//        for(int j = 0; j < D; j++)
//            KE += velocity(i,j) * velocity(i,j);
//
//        KE = 0.5 * KE;
//        cout<<"KE = "<<KE <<"\t PE = "<<PE<<"\t Total Energy = "<<PE + KE<<endl;
//        cout<<endl;
//        cout<<velocity;
//        cout<<endl;

        t +=dt;
        a++;
    }
*/


//    position_old(1,0) = 0.5;
//    position_old(1,1) = L/2;
//    position_old(2,0) = L-0.5;
//    position_old(2,1) = L/2;
//    position_old(3,0) = L/2;
//    position_old(3,1) = 0.5;
//    position_old(4,0) = L/2;
//    position_old(4,1) = L-0.5;
//
//    velocity(1,0) = -0.5;
//    velocity(1,1) = 0;
//    velocity(2,0) = 0.5;
//    velocity(2,1) = 0;
//    velocity(3,0) = 0;
//    velocity(3,1) = -0.5;
//    velocity(4,0) = 0;
//    velocity(4,1) = -0.5;
//    for(int j = 0; j < D; j++)
//        position_old(0,j) = L/2;

//    for(int i = 1; i < N; i++)
//        for(int j = 0; j < D; j++)
//        {
//            position_old(i, 0) = position_old(0,0) + i*cos(i*0.25 * PI);
//            position_old(i, 1) = position_old(0,1) + i*sin(i*0.25 * PI);
//        }

/*
    get_acceleration_2(position_old, acceleration, dist, PE, L);
    velocity = velocity + dt * acceleration;
    position = position_old + dt * velocity + 0.5 * dtsq * acceleration;
    periodic_boundary_condition(position, L);

    cout<<"position_old "<<endl<<position_old;
    cout<<"position "<<endl<<position;
    cout<<"velocity "<<endl<<velocity;
    cout<<"acceleration "<<endl<<acceleration;
    cout<<"distances "<<endl<<dist;
    cout<<"The minimum distance is "<<min_dist(dist);
    getch();
*/
