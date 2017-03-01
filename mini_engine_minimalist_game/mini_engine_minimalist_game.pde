// "mini" demo for minimalist game design class on 2/1/2017
// andy nealen, NYU, 2/1/2017
//
// delete the dots! 
// (with physics!) 
//
// simple euler method/time integration demo with 
// particles /w springs + poisson point sampling 
//
// CONTROLS/KEYS
// hit space bar to freeze time. hit it again to delete dots in the red ball 
// and continue the physics simulation
// default game is two tries of this (see maxDeletions), then 
// lowest remaining number of dots wins
//
// mr to reroll/reset the point sampling
// v to toggle velocity vector rendering
// s to toggle spring rendering
// p to pause/unpause game
// P to toggle poisson sampling (on by default)
// q to quit

// some GAME SPACE (parameters)
// max uses of the freeze/delete mechanic per playthrough/frame
int maxDeletions = 3; 
// radius of point freeze/deletion
float freezeRadius = 100.0;
boolean drawFreezeRadius = true;
// num points 
int NP = 101;
// spring constant/stiffness
float springStiffness = 0.25;
// min dist for particle-particle interaction
// for attractive forces
float D  = 150.0f;
// for repulsive forces
float rd = 60.0f;
// player particle mass
float pMass = 1.0;
// num frames (as in: bowling)
int numFrames = 5;
int frameNumber = 0; // runs from 0 to numFrames-1
int[] scores = new int[numFrames];
boolean renderScoreCard = true;

// IDEA: make a game with a single spring as a hand that grabs, and let it 
//       be flung by the main/heavy body (())----() -> (())-() -> ()---(())
//       skill lies in figuring out when to let go. could be a multiplayer
//       1-button 1D race game.

// IDEA: add in an adversarial NPC character who messes with the world in 
//       unexpected ways if "triggered": changes attractors, becomes a repulsor, etc.
//       thus introducing an obstacle that can/should not be deleted
//       adds some balance to the possibly overpowered delete mechanic
//       could also introduce a game end condition: weaken the adversaries 
//       links so it changes from indestructible to edible 
//       (yeah, i'm making ****ing osmos 2 :)

// instead of using a particle and particle system class
// just use vectors where index i refers to the same particle
// this should probably be implemented cleaner in production codeq
// but doesn't matter during prototyping!
PVector[] p = new PVector[NP]; // positions
PVector[] v = new PVector[NP]; // velocities
PVector[] f = new PVector[NP]; // forces (per timestep)
boolean[] alive = new boolean[NP]; // alive/dead status per particle

// global game states (again, somewhat bad practice... during prototyping it's ok!)
enum gs {
  GS_RUNNING, GS_FREEZE, GS_END
};
gs gameState; // game state
// used for player interaction -- only advance game state if invoked
boolean advanceState = false;
int numDeletions; // keep track of number of times freeze/delete was used

// boolean states for various render/sim options
boolean gameRunning = false; // for (un)pausing the simulation
boolean poissonSampling = true; // optional poisson sampling
boolean renderVelocities = false; // debug output
boolean renderSprings = true; // debug output
boolean renderFPS = true; // debug output
// last player positions
int PT = 20;
PVector[] trail = new PVector[PT]; 
int trailHead = 0; // index into trail[] to indicate the head of the trail (2d pos)
int drawTrail = 2; // (0 = no trail, 1 = not tapered, 2 = tapered trail
// last enemy positions
PVector[] eTrail = new PVector[PT]; 
boolean enemyPissedOff = false;
// particle rendering
boolean renderParticles = true;

// background color
color bg = color(0, 0, 0);

// random seed
int seed = 0;

// timer stuff
long timeStarted;
long timeElapsed;
long timeLastRenderTime;
long timeLastUpdateFPSTime;
long timeUpdateFPS = 1000; // update FPS counter every 1000 ms
long fps = 60;
long fpsSum = 0;
int numFPSSamples = 0; // for averaging the FPS over timeUpdateFPS
long timeTargetFrameTime = 8; // target render frametime in ms
long timeLastPhysicsTime;
PFont font;
boolean initPhysics = true; // start physics engine with a dt of 0, see below
long freezeTimer;
long freezeTime = 400;
boolean freezeTimerRunning = false;
float screenShakeIntensity = 0.08 ;
boolean useFreezeTimer = true;

// implemented a way to cycle through 4 speeds of the game
// without compromising on the physics stability
int numPhysicsSteps = 0;
int maxPhysicsSteps = 4;
float friction;

void setup() { // this gets called once

  // some stock setup code (screen, colors, linewidth
  size(600, 600, P2D);
  smooth();
  fill(0);
  stroke(255, 255, 255, 255);
  background(bg);
  strokeWeight(1);

  // timers
  font = createFont("Arial", 12, true);
  timeStarted = System.currentTimeMillis();
  timeLastRenderTime = timeStarted;
  timeLastUpdateFPSTime = timeStarted;
  timeLastPhysicsTime = timeStarted;

  // set up initial playfield w points
  resetGame();
} 

void draw() { // this gets called every frame of animation
  // use to also update the game state, not just render it
  // most game engines will have separate routines for
  // updating and rendering game state

  // if game is running and advanceState was invoked, update game logic
  // (see keyPressed() to see how advanceState is set to true)
  if (gameRunning && advanceState) {
    gameLogic();
    advanceState = false;
  }

  //freezeTimerRunning = false;
  // if the freeze state is induced, check to see if freezeTime has passed
  if (freezeTimerRunning) {
    long ct = System.currentTimeMillis();
    // check to see if freezeTimer has run out
    if (ct - freezeTimer > freezeTime) {
      // if yes, advance the game (to successor state) and stop the freeze timer
      gameLogic();
      freezeTimerRunning = false;
    }
  }

  // (physics) move and collide stuff if sim is running and is is the running state
  // this can always run at the highest framerate (= precision), but must be adaptive to the
  // framerate to be independent (minus numerical errors) of the rate at 
  // which draw() is running
  long deltaTime = System.currentTimeMillis() - timeLastPhysicsTime;
  if (gameRunning && gameState == gs.GS_RUNNING) {
    for (int m = 0; m < numPhysicsSteps + 1; m++) {
      move(deltaTime);
      collide();
    }
  }
  timeLastPhysicsTime = System.currentTimeMillis();

  // draw scene (all global vars for this prototype) even if not running game/logic
  // note: this is often locked to the refesh rate of the monitor in game engines
  // update 2/6/2017: added only to render if needed (see timeTargetFrameTime)
  if (renderTickPassed()) {
    renderStuff();
    long ct = System.currentTimeMillis();
    // compute current fps
    long dt = ct - timeLastRenderTime;
    fpsSum += 1000/dt;
    numFPSSamples++;
    // update global FPS once every timeUpdateFPS (to avoid display jitter)
    if (updateFPSTickPassed()) {
      fps = fpsSum / numFPSSamples;
      timeLastUpdateFPSTime = ct;
      numFPSSamples = 0;
      fpsSum = 0;
    }
    // output render framerate (if on)
    if (renderFPS) {
      textFont(font, 12);
      fill(255);
      text("FPS: " + fps, 10, 20);
      fill(0);
    }
    // update last render time
    timeLastRenderTime = ct;

    textFont(font, 12);
    fill(255);
    text("freezes remaining: " + numDeletions, 100, 20);
    text("particles remaining: " + numParticlesAlive(), 270, 20);
    text("speed: " + (numPhysicsSteps + 1) + "x", 450, 35);
    String F = String.format("%1$.6f", friction);
    text("friction: " + F, 450, 20);
    fill(0);
  }
}

// used to lock the framerate if so desired 
boolean renderTickPassed() {
  return (System.currentTimeMillis() - timeLastRenderTime) >= timeTargetFrameTime;
}

boolean updateFPSTickPassed() {
  return (System.currentTimeMillis() - timeLastUpdateFPSTime) >= timeUpdateFPS;
}

void resetGame() {
  genPoints();
  numDeletions = maxDeletions;
  gameState = gs.GS_RUNNING;
  println("GAME_STATE: RUNNING");
  println("remaining tries = " + numDeletions);  
  // init player and enemy trail
  for (int i = 0; i < PT; i++) {
    trail[i] = new PVector(p[0].x, p[0].y);
    eTrail[i] = new PVector(p[1].x, p[1].y);
  }
  enemyPissedOff = false;
  numPhysicsSteps = 0;
}

void resetGameAndScore() {
  resetGame();
  frameNumber = 0;
}

// render game based on debug flags and game state
void renderStuff() {

  // turn off the freeze timer if it's not being used (hack!)
  if (!useFreezeTimer) freezeTimerRunning = false;
  
  // screenshake during freeze state
  if (freezeTimerRunning) {
    pushMatrix();
    float shake = 100.0 * screenShakeIntensity;
    translate(random(-shake, shake), random(-shake, shake));
  }
  
  // clear framebuffer (background)
  background(bg);

  // draw connectors (if both particles are alive)
  if (renderSprings) {
    for (int i = 0; i < NP; i++) {
      for (int j = 0; j < NP; j++) {
        if (alive[i] && alive[j]) {
          float d = p[i].dist(p[j]);
          if (d < D) {
            // set stroke alpha (transparency) based on distance threshold
            if (enemyPissedOff && 
               (p[1].dist(p[i]) < freezeRadius || p[1].dist(p[j]) < freezeRadius) &&
                i != 1 && j != 1) {
              stroke(255, 0, 0, D-d);
            } else {
              stroke(255, 255, 255, D-d);
            }
            line(p[i].x, p[i].y, p[j].x, p[j].y);
            // reset alpha to 1/255
            stroke(255, 255, 255, 255);
          }
        }
      }
    }
  }

  // draw velocity vectors (if alive) 
  if (renderVelocities) {
    float velRenderScale = 1.0;
    for (int i = 0; i < NP; i++) {
      if (alive[i]) {
        stroke(255, 255, 255, 255);
        line(p[i].x, p[i].y, 
          p[i].x + velRenderScale*v[i].x, 
          p[i].y + velRenderScale*v[i].y);
        stroke(255, 255, 255, 255);
      }
    }
  }

  // draw points (if alive)
  if (renderParticles) {
    for (int i = 1; i < NP; i++) {
      if (alive[i]) {
        ellipse(p[i].x, p[i].y, 10, 10);
      }
    }
  }

  // draw player and enemy particle with or without trail 
  if (drawTrail < 1) {
    fill(255, 255, 255, 255);
    stroke(255, 255, 255, 255);
    ellipse(p[0].x, p[0].y, 20, 20);
    ellipse(p[1].x, p[1].y, 20, 20);
    fill(150, 0, 0, 255);
    stroke(150, 0, 0, 255);
    ellipse(p[0].x, p[0].y, 15, 15);
    if (!enemyPissedOff) {
      fill(0, 0, 150, 255);
      stroke(0, 0, 150, 255);
    }
    ellipse(p[1].x, p[1].y, 15, 15);
  } else {
    int size = 20;
    int trailIndex = trailHead;
    fill(255, 255, 255, 255);
    stroke(255, 255, 255, 255);
    for (int i = 0; i < PT; i++) {
      trailIndex--; 
      if (trailIndex < 0) trailIndex = PT - 1;
      PVector pos = trail[trailIndex];
      ellipse(pos.x, pos.y, size, size);
      pos = eTrail[trailIndex];
      ellipse(pos.x, pos.y, size, size);
      if (drawTrail > 1) size *= 0.97;
    }
    size = 15;
    trailIndex = trailHead;
    fill(150, 0, 0, 255);
    stroke(150, 0, 0, 255);
    for (int i = 0; i < PT; i++) {
      trailIndex--; 
      if (trailIndex < 0) trailIndex = PT - 1;
      PVector pos = trail[trailIndex];
      ellipse(pos.x, pos.y, size, size);
      if (drawTrail > 1) size *=0.94;
    }
    size = 15;
    trailIndex = trailHead;
    if (!enemyPissedOff) {
      fill(0, 0, 150, 255);
      stroke(0, 0, 150, 255);
    }
    for (int i = 0; i < PT; i++) {
      trailIndex--; 
      if (trailIndex < 0) trailIndex = PT - 1;
      PVector pos = eTrail[trailIndex];
      ellipse(pos.x, pos.y, size, size);
      if (drawTrail > 1) size *=0.94;
    }
  }

  fill(0);
  stroke(255, 255, 255, 255);

  // if game is paused, render alpha blended quad over scene
  if (!gameRunning) {
    fill(0, 0, 0, 150);
    rect(1, 1, 597, 597);
    fill(0);
  }

  // if gameState is FREEZE, show circle of deletion for better understanding
  //  if (gameState == gs.GS_FREEZE) {
  if (drawFreezeRadius) {
    fill(150, 0, 0, 100);
    stroke(150, 0, 0, 100);
    ellipse(p[0].x, p[0].y, freezeRadius*2, freezeRadius*2);
    fill(0);
    stroke(255, 255, 255, 255);
  }
  //  }

  // screenshake during freeze state (restore identity transformation with popMatrix())
  if (freezeTimerRunning) {
    popMatrix();
  } // so, no shake on anything after this line...

  if (renderScoreCard) renderScoreCard(10, 540, numFrames);
}

void renderScoreCard(int x, int y, int numScores) {
  fill(255);
  textFont(font, 12);
  text("SCORECARD", x, y);
  // horizontal lines
  line(x, y + 10, x + 30 * (numScores +1), y + 10);
  line(x, y + 30, x + 30 * (numScores +1), y + 30);
  line(x, y + 50, x + 30 * (numScores +1), y + 50);
  // vertical lines and frame labels/numbers
  line(x, y + 10, x, y + 50);
  int i;
  for (i = 0; i < numScores; i++) {
    text(i + 1, x + 13 + (i * 30), y + 25);
    line(x + 30 * (i + 1), y + 10, x + 30 * (i + 1), y + 50);
    if (frameNumber > i) {
      text(scores[i], x + 13 + (i * 30), y + 45);
    }
  }
  text("T", x + x + (i * 30), y + 25);
  if (frameNumber == numFrames) {
    int score = 0;
    for (int s = 0; s < numFrames; s++) {
      score += scores[s];
    }
    text(score, x + x + (i * 30), y + 45);
  } 
  line(x + 30 * (i + 1), y + 10, x + 30 * (i + 1), y + 50);
}

void keyPressed() {

  // don't react to key presses while freeze timer is running
  if (freezeTimerRunning) return;

  if (key == ' ') {
    // notify game loop that player interaction occurred
    // (see gameLogic() in draw() above)
    if (gameRunning) {
      advanceState = true;
    } else {
      gameRunning = true;
    }
  } else if (key == 'p') {
    // (un)pause
    gameRunning = !gameRunning;
  } else if (key == 'r') {
    seed = 0;
    resetGameAndScore();
  } else if (key == 'q') {
    exit();
  } else if (key == 'P') {
    println();
    poissonSampling = !poissonSampling; 
    if (poissonSampling) {
      println("POISSON SAMPLING");
    } else {
      println("RANDOM SAMPLING");
    }
  } else if (key == 'v') {
    renderVelocities = !renderVelocities;
  } else if (key == 's') {
    renderSprings = !renderSprings;
  } else if (key == 'S') {
    renderScoreCard = !renderScoreCard;
  } else if (key == 'F') {
    useFreezeTimer = !useFreezeTimer;
  } else if (key == 'a') {
    if (screenShakeIntensity < 0.01) {
      screenShakeIntensity = 0.08;
    } else if (screenShakeIntensity > 0.05) {
      screenShakeIntensity = 0.0;
    }
  } else if (key == 'f') {
    renderFPS= !renderFPS;
  } else if (key == 't') {
    drawTrail += 1;
    drawTrail %= 3;
  } else if (key == 'y') {
    renderParticles = !renderParticles;
  } else if (key == 'e') {
    drawFreezeRadius = !drawFreezeRadius;
  } else if (key == 'b') {
    numPhysicsSteps++;
    numPhysicsSteps %= maxPhysicsSteps;
  }
}

// this is the "big" switch/case statement for all game state based logic
// this holds the logic "flow" of the game. pretty "hacky," and can be improved
// (does not matter for small prototypes)
// cycles through RUNNING -> FREEZE -> RUNNING
// unless the game needs to be reset (all deletes used) and then runs
// RUNNING -> FREEZE -> END -> RUNNING
// in order to compute and display score, and then reset
void gameLogic() {
  switch (gameState) {
  case GS_RUNNING:
    gameState = gs.GS_FREEZE;
    // after setting game state to freeze, start freeze timer
    freezeTimerRunning = true;
    freezeTimer = System.currentTimeMillis();
    println("GAME_STATE: FREEZE");
    break;
  case GS_FREEZE:
    // set all particles in range to no longer alive
    for (int i = 2; i < NP; i++ ) {
      if (p[0].dist(p[i]) < freezeRadius) {
        alive[i] = false;
      }
    }
    // piss off/pacify enemy in freeze radius
    if (p[0].dist(p[1]) < freezeRadius) {
      enemyPissedOff = !enemyPissedOff;
    }
    // used up one deletion (see maxDeletions)
    numDeletions--;
    // state transition at the end of GS_FREEZE: running or end
    if (numDeletions > 0 && numParticlesAlive() > 0) {
      gameState = gs.GS_RUNNING;
      println("GAME_STATE: RUNNING");
      println("remaining tries = " + numDeletions);
    } else {
      gameState = gs.GS_END;
      println("GAME_STATE: END");
      println("particles remaining: " + numParticlesAlive());
      if (frameNumber < numFrames) {
        // store current score
        scores[frameNumber] = numParticlesAlive();
        frameNumber++;
      }
    }
    break;
  case GS_END:
    if (frameNumber == numFrames) frameNumber = 0;
    resetGame();
    break;
  }
}

// all physics motion without collision correction for walls (see collision() below)
void move(long deltaTime) { // deltaTime in milliseconds

  // so the engine doesn't do weird stuff on frame one due to too large a dt
  if (initPhysics) {
    deltaTime = 0;
    // only do this once
    initPhysics = false;
  }

  // clear (temporary, per frame of animation) forces
  for (int i = 0; i < NP; i++) {
    f[i] = new PVector(0, 0);
  }

  // add all particle-particle interaction forces
  for (int i = 0; i < NP; i++) {
    for (int j = 0; j < NP; j++) {
      if (alive[i] && alive[j]) {
        PVector f_ij = PVector.sub(p[j], p[i]); // f_ij = vecotr with length dist(i,j)

        float dist = p[i].dist(p[j]);
        // if these two particles are close...
        if (dist < D) {

          // feedback loop based on fraction of particles left?
          // this effectively initially decreases, then increases the spring stiffness
          float fractionAlive = (float) numParticlesAlive() / (float) NP;
          float adjustedStiffness = springStiffness * (0.6f + 0.4 * (1.0f - fractionAlive));
          // or don't do this
          // float adjustedStiffness = springStiffness;

          // ... pull a tiny bit towards other particle, proportional to
          // springStiffness times particle distance (= length of "spring" multiplied
          // by the springStiffness constant
          f_ij.mult(adjustedStiffness); // f_ij = springStiffness * dist(i,j)
          f[i].add(f_ij);
          
          // if the enemy is pissed off, go insane (i.e. mess with other particle forces) 
          if (enemyPissedOff && 
             (p[1].dist(p[i]) < freezeRadius || p[1].dist(p[j]) < freezeRadius) &&
              i != 1 && j != 1) {
            f_ij.mult(2.0);
            f[i].sub(f_ij);
          }

          // some repulsion if particles are too close
          if (dist < rd) {
            PVector r_ij = PVector.sub(p[j], p[i]); // r_ij
            r_ij.normalize();
            r_ij.mult((rd - dist) * springStiffness);
            f[i].sub(r_ij);
          }

          // super weird rubber effect instead 
          // (repulsion, makes game worse, looks cool. use instead of fp[i].add(f_ij) above)
          // f[i].sub(f_ij);
        }
      }
    }
  }

  // euler method per particle: update v[i] and p[i] based on f[i]
  float dt = deltaTime/1000.0f; // delta time (dt) in seconds
  float m  = 1.0; // mass of particles in kg
  // newton's second law of motion: f(orce) = m(ass) * a(cceleration, or a = f / m
  // assuming a timestep dt (in seconds), the difference in velocity 
  // is dv = a * dt, so the updated velocity is v(next)= v(current) + dv
  // the difference in position is dp = v * dt and so p(next) = p(current) + dp = p + v * dt
  //
  // we're "simulating" a constant mass (1.0) but not timestep
  // most engines let the "update" function know how much time has passed
  // such that dt can be variable, and the physics can adjust to how
  // much real-time has passed, as is the case here
  // note that this is a bit overkill for a prototype, but needed for production

  // add in some (game state dependent) friction
  float fractionAlive = (float) (numParticlesAlive() + 2) / (float) NP;
  friction = 0.996f + 0.004f * fractionAlive;
  // or not
  // friction = 1.0f;

  // perform one (explicit) euler timestep
  for (int i = 0; i < NP; i++) {
    m = 1.0; 
    if (i == 0) m = pMass;   // mass is 1kg by default, player mass otherwise
    f[i].div(m);                      // a[i]   = f[i] / m[i]      (f  = m * a)
    f[i].mult(dt);                    // d_v[i] = a[i] * dt        (dv = a * dt)
    v[i].add(f[i]);                   // v[i]   = v[i] + d_v[i]    (v  = v + dv)
    v[i].mult(friction);              // v[i]   = v[i] * f         (friction dampens velocity)
    p[i].add(PVector.mult(v[i], dt));  // p[i]   = p[i] + v[i] * dt (p  = p + v * dt)
  }

  // store player (physics) positions in an array for trail rendering
  trail[trailHead] = new PVector(p[0].x, p[0].y);
  eTrail[trailHead] = new PVector(p[1].x, p[1].y);
  trailHead++;
  trailHead = trailHead % PT;
}

// contains all collision corrections, here only with the wall, but could be extended
// to contain and react to particle-particle collisions/overlaps
void collide() {

  // check wall collisions  
  for (int i = 0; i < NP; i++) {
    if (alive[i]) {
      if (p[i].x > width) {
        v[i].x = -v[i].x;
        p[i].x = width;
      }
      if (p[i].x < 0) {
        v[i].x = -v[i].x;
        p[i].x = 0;
      }
      if (p[i].y > height) {
        v[i].y = -v[i].y;
        p[i].y = height;
      }
      if (p[i].y > height || p[i].y < 0) {
        v[i].y = -v[i].y;
        p[i].y = 0;
      }
    }
  }
}

// score: how many particles are alive?
int numParticlesAlive() {
  int s = 0;
  for (int i = 2; i < NP; i++) {
    if (alive[i]) {
      s++;
    }
  }
  return s;
}

// everything below here are routines to generate and distribute points in the plane

PVector genPoint() {
  return new PVector(random(20, width -20), random(20, height-20));
}

void genPoints() {
  randomSeed(seed);
  if (poissonSampling) {
    rerollPoisson();
  } else {
    reroll();
  }
  seed++;
}

void reroll() {
  for (int i = 0; i < NP; i++) {
    p[i] = genPoint();
    // set initial velocity to zero
    v[i] = new PVector(0, 0);
    // initially alive
    alive[i] = true;
  }
}

void rerollPoisson() { // dart throwing algorithm with... 

  // ... maxTries and minDist
  int maxTries = 100;
  float minDist = 40.0;

  println("------------------------------");
  println("new poisson run with maxTries = " + maxTries);

  // no tries yet...
  int numTries = 0;

  for (int i = 0; i < NP; ) {

    // grab a random point
    PVector pt = genPoint();
    // assume it's valid
    boolean validPoint = true;
    // increment attempt num
    numTries++;

    // test for validity of point
    for (int j = 0; j < i; j++) {
      float dist = PVector.dist(pt, p[j]); 
      if (dist < minDist) {
        // too close to existing point, reject sample point
        validPoint = false;
        // we can stop checking. a single fail is a fail.
        break;
      }
    }

    // if point is valid, or we maxed out numTries, use it 
    if (validPoint || numTries > maxTries) {
      p[i] = pt;
      // set initial velocity to zero
      v[i] = new PVector(0, 0);
      // initially alive
      alive[i] = true;
      // if it took more than one attempt, output how many it took
      if (numTries > 1) {
        println("point " + i + " rejected, num tries to find: " + numTries);
      }
      // reset: go to next point sample, reset try counter
      i++;
      numTries = 0;
    }
  }  
  println();
}