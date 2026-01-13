class Simulation{

    constructor(){
        this.particles = [];
        this.particleEmitters = [];
        this.springs = new Map();
        this.shapes = [];

        this.numParticles = 1000;
        this.velocityDamping = 0.99;
        this.gravity = new Vector2(0, 1);
        this.restDensity = 10;
        this.kNear = 2;
        this.k = 0.3;
        this.interactionRadius = 20;

        //viscosity
        this.sigma = 0.1;
        this.beta = 0;

        //plasticity
        this.gamma = 0.3;
        this.plasticity = 0.7;
        this.springStiffness = 0.3;

        this.fluidHashGrid = new FluidHashGrid(this.interactionRadius);
        this.instantiateParticles();
        this.fluidHashGrid.initialize(this.particles);

        this.emitter = this.createParticleEmitter(
            new Vector2(canvas.width / 2, 400), // position
            new Vector2(0,-1), // direction
            20, // size
            1.5,  // spawn interval (increased)
            5, // amount (reduced)
            20  // speed
        );
        
        this.emitterEnabled = false;


    }

    createParticleEmitter(position, direction, size, spawnInterval, amount, velocity){
        let emitter = new ParticleEmitter(position, direction, size, spawnInterval, amount, velocity);
        this.particleEmitters.push(emitter);
        return emitter;
    }

    getShapeAt(pos){
        // Check if clicking on emitter
        let distToEmitter = Sub(pos, this.emitter.position).Length();
        if(distToEmitter < 15){
            return this.emitter;
        }
        
        // Check shapes
        for(let i=0; i< this.shapes.length; i++){
            if(this.shapes[i].isPointInside(pos)){
                return this.shapes[i];
            }
        }
        return null;
    }

    instantiateParticles(){
        let offsetBetweenParticles = 10;
        let offsetAllParticles = new Vector2(canvas.width / 2 - 400, canvas.height / 2 - 400);

        let xParticles = Math.sqrt(this.numParticles);
        let yParticles = xParticles;

        for(let x=0; x< xParticles; x++){
            for(let y=0; y< yParticles; y++){
                let position = new Vector2(
                    x*offsetBetweenParticles +offsetAllParticles.x,
                    y*offsetBetweenParticles +offsetAllParticles.y);
                
                let particle = new Particle(position);
                //particle.velocity = Scale(new Vector2(-0.5 + Math.random(), -0.5 + Math.random()), 200);
                this.particles.push(particle);
            }

        }
    }

    neighbourSearch(){
        this.fluidHashGrid.clearGrid();
        this.fluidHashGrid.mapParticleToCell();
    }

    update(dt){
        this.neighbourSearch();
        if(this.emitterEnabled){
            this.emitter.spawn(dt, this.particles);
        }
        if(this.rotate){
            this.emitter.rotate(0.01);
        }
        this.applyGravity(dt);
        this.viscosity(dt);
        this.predictPositions(dt);
        this.adjustSprings(dt);
        this.springDisplacement(dt);
        this.doubleDensityRelaxation(dt);
        this.worldBoundary();
        this.handleShapeCollisions(dt);
        this.computeNextVelocity(dt);


    }

    adjustSprings(dt){
        for(let i=0; i< this.particles.length; i++){
            let neighbourIndices = this.fluidHashGrid.getNeighbourOfParticleIdx(i);
            let particleA = this.particles[i];

            for(let j = 0; j < neighbourIndices.length;j++){
                let neighbourIdx = neighbourIndices[j];
                if(i == neighbourIdx) continue;

                let particleB = this.particles[neighbourIdx];
                let springId = i < neighbourIdx ? i + neighbourIdx * this.particles.length : neighbourIdx + i * this.particles.length;

                if(this.springs.has(springId)){
                    continue;
                }

                let rij = Sub(particleB.position,particleA.position); 
                let q = rij.Length() / this.interactionRadius;

                if(q < 1){
                    let newSpring = new Spring(i, neighbourIdx, this.interactionRadius);
                    this.springs.set(springId, newSpring);
                }
            }
        }


        for(let [key, spring] of this.springs){
            let pi = this.particles[spring.particleAIdx];
            let pj = this.particles[spring.particleBIdx];

            let rij = Sub(pi.position, pj.position).Length();
            let Lij = spring.length;
            let d = this.gamma * Lij;

            if(rij > Lij + d){
                spring.length += dt * this.plasticity * (rij - Lij - d); // stretching

            }else if(rij < Lij - d){ 
                spring.length -= dt * this.plasticity * (Lij - d - rij); // compression
            }

            if(spring.length > this.interactionRadius){
                this.springs.delete(key);
            }
        }
    }

    springDisplacement(dt){
        let dtSquared = dt * dt;

        for(let [key, spring] of this.springs){
            let pi = this.particles[spring.particleAIdx];
            let pj = this.particles[spring.particleBIdx];

            let rij = Sub(pi.position, pj.position);
            let distance = rij.Length();

            if(distance < 0.0001){
                continue;
            }

            rij.Normalize();
            let displacementTerm = dtSquared * this.springStiffness * 
                (1 - spring.length / this.interactionRadius) * (spring.length - distance);

            rij = Scale(rij, displacementTerm * 0.5);

            pi.position = Add(pi.position, rij);
            pj.position = Sub(pj.position, rij);
        }
    }


    viscosity(dt){
        for( let i=0; i< this.particles.length; i++){
            let neighbourIndices = this.fluidHashGrid.getNeighbourOfParticleIdx(i);
            let particleA = this.particles[i];

            for( let j=0; j< neighbourIndices.length; j++){
                let particleB = this.particles[neighbourIndices[j]];
                if(i == neighbourIndices[j]){
                    continue;
                }
                let rij = Sub(particleB.position, particleA.position);
                let velocityA = particleA.velocity;
                let velocityB = particleB.velocity;
                let q = rij.Length() / this.interactionRadius;

                if(q < 1){
                    rij.Normalize();
                    let u = Sub(velocityA, velocityB).Dot(rij);

                    if(u > 0){
                        let ITerm = dt * (1 - q) * (this.sigma * u + this.beta * u * u);
                        let I = Scale(rij, ITerm);

                        particleA.velocity = Sub(particleA.velocity, Scale(I, 0.5));
                        particleB.velocity = Add(particleB.velocity, Scale(I, 0.5));
                    }
                }
            }
        }
    }

    doubleDensityRelaxation(dt){
        for( let i=0; i< this.particles.length; i++){
            let density = 0;
            let densityNear = 0;
            let neighbourIndices = this.fluidHashGrid.getNeighbourOfParticleIdx(i);
            let particleA = this.particles[i];

            for( let j=0; j< neighbourIndices.length; j++){
                let particleB = this.particles[neighbourIndices[j]];
                if(i == neighbourIndices[j]){
                    continue;
                }
                let rij = Sub(particleB.position, particleA.position);
                let q = rij.Length() / this.interactionRadius;

                if(q < 1){
                    density += Math.pow(1-q, 2);
                    densityNear += Math.pow(1-q, 3);
                }
            }

            let pressure = this.k * (density - this.restDensity);
            let pressureNear = this.kNear * densityNear;
            let particleADisplacement = Vector2.Zero();

            for( let j=0; j< neighbourIndices.length; j++){
                let particleB = this.particles[neighbourIndices[j]];
                if(i == neighbourIndices[j]){
                    continue;
                }
                let rij = Sub(particleB.position, particleA.position);
                let q = rij.Length() / this.interactionRadius;

                if(q < 1){
                    rij.Normalize();
                    let displacementTerm = Math.pow(dt, 2) * (pressureNear * (1-q) + Math.pow(1-q, 2)); 
                    let D = Scale(rij, displacementTerm);
                    particleB.position = Add(particleB.position, Scale(D, 0.5));
                    particleADisplacement = Sub(particleADisplacement, Scale(D, 0.5));
                }
            }
            particleA.position = Add(particleA.position, particleADisplacement);
        }
    }

    applyGravity(dt){
        for( let i=0; i< this.particles.length; i++){
            this.particles[i].velocity = Add(this.particles[i].velocity, Scale(this.gravity, dt));
        }
    }

    predictPositions(dt){
        for(let i=0; i< this.particles.length; i++){
            this.particles[i].prevPosition = this.particles[i].position.Cpy();
            let positionDelta = Scale(this.particles[i].velocity, dt * this.velocityDamping);
            this.particles[i].position = Add(this.particles[i].position, positionDelta);
        }
    }

    computeNextVelocity(dt){
        for(let i=0; i< this.particles.length; i++){
            let velocity = Scale(Sub(this.particles[i].position, this.particles[i].prevPosition), 1.0 /dt);
            this.particles[i].velocity = velocity;

        }
    }

    worldBoundary(){
        for(let i=0; i< this.particles.length; i++){
            let pos = this.particles[i].position;

            if(pos.x < 0){
                this.particles[i].position.x = 0;
                this.particles[i].prevPosition.x = 0;
            }
            if(pos.y < 0){
                this.particles[i].position.y = 0;
                this.particles[i].prevPosition.y = 0;
            }
            if(pos.x > canvas.width){
                this.particles[i].position.x = canvas.width-1;
                this.particles[i].prevPosition.x = canvas.width-1;
            }
            if(pos.y > canvas.height){
                this.particles[i].position.y = canvas.height-1;
                this.particles[i].prevPosition.y = canvas.height-1;
            }
        }
    }

    handleShapeCollisions(dt){
        for(let i = 0; i < this.particles.length; i++){
            let particle = this.particles[i];
            
            for(let s = 0; s < this.shapes.length; s++){
                let shape = this.shapes[s];
                
                if(shape.isPointInside(particle.position)){
                    let collisionData = this.getCollisionData(particle.position, shape);
                    
                    if(collisionData){
                        let normal = collisionData.normal;
                        let penetration = collisionData.penetration;
                        
                        // Push particle out of shape
                        particle.position = Add(particle.position, Scale(normal, penetration + 0.1));
                        
                        // Dampen velocity and reflect component along normal
                        let vel = Sub(particle.position, particle.prevPosition);
                        let velAlongNormal = vel.Dot(normal);
                        
                        if(velAlongNormal < 0){
                            // Remove velocity component along normal and dampen it
                            let normalComponent = Scale(normal, velAlongNormal);
                            vel = Sub(vel, normalComponent);
                            vel = Sub(vel, Scale(normalComponent, 0.8)); // friction/damping
                            particle.prevPosition = Sub(particle.position, vel);
                        }
                    }
                }
            }
        }
    }
    
    getCollisionData(particlePos, shape){
        if(shape instanceof Circle){
            let toParticle = Sub(particlePos, shape.position);
            let distance = toParticle.Length();
            
            if(distance < 0.001) distance = 0.001; // avoid division by zero
            
            let normal = Scale(toParticle, 1 / distance);
            let penetration = shape.radius - distance;
            
            return { normal: normal, penetration: penetration };
        }
        else if(shape instanceof Polygon){
            return this.getPolygonCollisionData(particlePos, shape);
        }
        
        return null;
    }
    
    getPolygonCollisionData(particlePos, polygon){
        let minDistance = Infinity;
        let closestEdgeIdx = -1;
        let closestPointOnEdge = null;
        
        // Check distance to each edge
        for(let i = 0; i < polygon.vertices.length; i++){
            let v1 = polygon.vertices[i];
            let v2 = polygon.vertices[(i + 1) % polygon.vertices.length];
            
            let edge = Sub(v2, v1);
            let toParticle = Sub(particlePos, v1);
            
            // Project particle onto edge
            let edgeLenSq = edge.Length2();
            let t = Math.max(0, Math.min(1, toParticle.Dot(edge) / edgeLenSq));
            
            let closest = Add(v1, Scale(edge, t));
            let dist = Sub(particlePos, closest).Length();
            
            if(dist < minDistance){
                minDistance = dist;
                closestEdgeIdx = i;
                closestPointOnEdge = closest;
            }
        }
        
        if(closestEdgeIdx >= 0 && polygon.normals && polygon.normals[closestEdgeIdx]){
            let normal = polygon.normals[closestEdgeIdx].Cpy();
            
            // Ensure normal points outward (away from polygon center)
            let toParticle = Sub(particlePos, polygon.vertices[closestEdgeIdx]);
            if(toParticle.Dot(normal) < 0){
                normal = Scale(normal, -1);
            }
            
            return { normal: normal, penetration: Math.max(0.5, 2 - minDistance) };
        }
        
        return null;
    }

    draw(){
        for(let i=0; i< this.shapes.length; i++){
            this.shapes[i].draw();
        }

        for(let i=0; i< this.particles.length; i++){
            let position = this.particles[i].position;
            let color = this.particles[i].color;
            DrawUtils.drawPoint(position, 5, color);
        }

        for(let i=0; i< this.particleEmitters.length; i++){
            this.particleEmitters[i].draw();
        }

    }

    toggleEmitter(){
        this.emitterEnabled = !this.emitterEnabled;
    }

    addCircle(x, y, radius){
        let circle = new Circle(new Vector2(x, y), radius, "orange");
        this.shapes.push(circle);
    }

    removeCircle(){
        for(let i = this.shapes.length - 1; i >= 0; i--){
            if(this.shapes[i] instanceof Circle){
                this.shapes.pop();
                break;
            }
        }
    }
}