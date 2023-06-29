const GM = 1.32712440018e2;	// Mm^3/kg/s^2
const AU =  1.495978707e5;	// Mm

var view = 3;
var init = 0;

const spacecraftsets = [
  {
    name: "Circular",
    tail: 300,
    states: [
      { x: AU-1, y: -1, z: -3, vx:1e-8, vy:30e-3, vz:2e-7, r:255, g:128, b:128 },
      { x: AU+1, y: -1, z: 2, vx:0, vy:30e-3-4e-7, vz:0, r:128, g:255, b:128 },
      { x: AU-1, y: 1, z: 5, vx:0, vy:30e-3, vz:2e-7, r:128, g:128, b:255 },
      { x: AU+1, y: 2, z: 2, vx:0, vy:30e-3-4e-7, vz:0, r:255, g:255, b:0 }
    ]
  },
  {
    name: "Eccentric",
    tail: 600,
    states: [
      { x: AU-1, y: -1, z: -1, vx:0e-7, vy:35e-3-15e-8, vz:1e-7, r:255, g:128, b:128 },
      { x: AU+1, y: -1, z: 1, vx:0e-7, vy:35e-3-50e-8, vz:0, r:128, g:255, b:128 },
      { x: AU-1, y: 1, z: -2, vx:0e-7, vy:35e-3-15e-8, vz:1e-7, r:128, g:128, b:255 },
      { x: AU+1, y: 1, z: 3, vx:0e-7, vy:35e-3-50e-8, vz:0e-7, r:255, g:255, b:0 }
    ]
  },
  {
    name: "High ecc.",
    tail: 600,
    states: [
      { x: 0.6*AU-0.5, y: -0.5, z: -0.5, vx: 1e-7, vy:0.0485+1.7e-7, vz: 5e-7, r:255, g:128, b:128 },
      { x: 0.6*AU+0.5, y: -0.5, z:  0.5, vx: 1e-7, vy:0.0485-1.7e-7, vz:-8e-7, r:128, g:255, b:128 },
      { x: 0.6*AU-0.5, y:  0.5, z:  0.5, vx: 1e-7, vy:0.0485+1.7e-7, vz:-4e-7, r:128, g:128, b:255 },
      { x: 0.6*AU+0.5, y:  0.5, z: -0.5, vx:-1e-7, vy:0.0485-1.7e-7, vz: 3e-7, r:255, g:255, b:0 }
    ]
  },
  {
    name: "High ecc. (tiny)",
    tail: 600,
    states: [
      { x: 0.6*AU-0.5e-1, y: -0.5e-1, z: -0.5e-1, vx: 1e-8, vy:0.0485+1.7e-8, vz: 5e-8, r:255, g:128, b:128 },
      { x: 0.6*AU+0.5e-1, y: -0.5e-1, z:  0.5e-1, vx: 1e-8, vy:0.0485-1.7e-8, vz:-8e-8, r:128, g:255, b:128 },
      { x: 0.6*AU-0.5e-1, y:  0.5e-1, z:  0.5e-1, vx: 1e-8, vy:0.0485+1.7e-8, vz:-4e-8, r:128, g:128, b:255 },
      { x: 0.6*AU+0.5e-1, y:  0.5e-1, z: -0.5e-1, vx:-1e-8, vy:0.0485-1.7e-8, vz: 3e-8, r:255, g:255, b:0 }
    ]
  },
  {
    name: "High ecc. (large)",
    tail: 600,
    states: [
      { x: 0.6*AU-0.5e1, y: -0.5e1, z: -0.5e1, vx: 1e-7, vy:0.0485+1.7e-7, vz: 5e-7, r:255, g:128, b:128 },
      { x: 0.6*AU+0.5e1, y: -0.5e1, z:  0.5e1, vx: 1e-7, vy:0.0485-1.7e-7, vz:-8e-7, r:128, g:255, b:128 },
      { x: 0.6*AU-0.5e1, y:  0.5e1, z:  0.5e1, vx: 1e-7, vy:0.0485+1.7e-7, vz:-4e-7, r:128, g:128, b:255 },
      { x: 0.6*AU+0.5e1, y:  0.5e1, z: -0.5e1, vx:-1e-7, vy:0.0485-1.7e-7, vz: 3e-7, r:255, g:255, b:0 }
    ]
  },
  {
    name: "test",
    tail: 300,
    states: [
      { x: AU, y: 0, z: 0, vx:0, vy:30e-3, vz:0, r:255, g:128, b:128 },
      { x: AU+2, y: 0, z: 0, vx:0, vy:30e-3, vz:4e-7, r:128, g:255, b:128 },
      { x: AU, y: 2, z: 0, vx:0, vy:30e-3, vz:0, r:128, g:128, b:255 },
      { x: AU, y: 0, z: 2, vx:0, vy:30e-3, vz:0, r:255, g:255, b:0 }
    ]
  },
];

class State
{
  constructor(x, y, z, vx, vy, vz, r, g, b)
  {
    this.x = x;
    this.y = y;
    this.z = z;
    this.vx = vx;
    this.vy = vy;
    this.vz = vz;
	this.r = r;
	this.g = g;
	this.b = b;
  }

  add(otherState)
  {
    let newX = this.x + otherState.x;
    let newY = this.y + otherState.y;
    let newZ = this.z + otherState.z;
    let newVx = this.vx + otherState.vx;
    let newVy = this.vy + otherState.vy;
    let newVz = this.vz + otherState.vz;

    return new State(newX, newY, newZ, newVx, newVy, newVz, this.r, this.g, this.b);
  }

  multiply(dt)
  {
    return new State(this.x*dt, this.y*dt, this.z*dt, this.vx*dt, this.vy*dt, this.vz*dt, this.r, this.g, this.b);
  }
}

function saveAll()
{
	let savedata = JSON.stringify(
	{
		init: init,
		view: view,
		states: states,
		time: time, top: 
		{
			size: topView.bufferSize,
			buffer:topView.buffer
		},
		rts: {M: rtsM, T: rtsT},
		camera: camera
	});

	let blob = new Blob([savedata], {type: "application/json"});
	let url = URL.createObjectURL(blob);
	let a = document.createElement("a");
	a.href = url;
	a.download=document.title + ".json";
	a.click();
}

function loadAll()
{
	stop();

	let fileInput = document.createElement("input");
	fileInput.type = "file";
	fileInput.accept = ".json, application/json";
	fileInput.addEventListener("change", () =>
	{
		let file = fileInput.files[0];
		fileName = file.name.replace(/.json$/,"");
		let reader = new FileReader();
		reader.addEventListener("load", () =>
		{
			document.getElementById("init").disabled = true;

			let data = JSON.parse(reader.result);
			init = data.init;
			view = data.view;
			time = data.time;
			topView.bufferSize = data.top.size;
			topView.buffer = data.top.buffer;
			rtsM = new TS(data.rts.M.size);
			rtsM.data = data.rts.M.data;
			rtsT = new TS(data.rts.T.size);
			rtsT.data = data.rts.T.data;
			camera.x = data.camera.x;
			camera.y = data.camera.y;
			camera.z = data.camera.z;

			document.getElementById("init").value = init;
			doResize();
			onInit(init);
			states = [];
			data.states.forEach(s =>
			{
				states.push(new State(s.x, s.y, s.z, s.vx, s.vy, s.vz, s.r, s.g, s.b));
			});
			document.getElementById("view").value = view;
			setView(view);
			fileInput.remove();
		});
		reader.readAsText(file);
	});
	fileInput.click();
}


const camera = { x: 0, y: 0, z: 10, phi: 0, theta: 0 };
var spacecrafts;
var states;
var indices;

// Force modifiers:
// w - induced rotation, can test with 1e-9
// m - Yukawa mass (inverse range) in 1/AU
// y - Yukawa coupling constant
var MOD = { w: 0, m: 0, y: 0 };

function ag(s)
{
	var r = norm(s); // Math.sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
	var r3 = r * r * r;
	var GMY = GM * (1 + MOD.y * (1 - (1 + r*MOD.m/AU)*Math.exp(-r*MOD.m/AU)));
	return new State(0, 0, 0, -GMY * s.x / r3, -GMY * s.y / r3, -GMY * s.z / r3);
}

function a(s)
{
  let a0 = ag(s);

  let ex = states.reduce((a,b)=>({x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}));
  let nx = norm(ex);
  ex.x /= nx;
  ex.y /= nx;
  ex.z /= nx;

  let r = norm(s);

  a0.vx += MOD.w/r * (ex.y*s.z - ex.z*s.y);
  a0.vy += MOD.w/r * (ex.z*s.x - ex.x*s.z);
  a0.vz += MOD.w/r * (ex.x*s.y - ex.y*s.x);

  return a0;
}

function rk4(s, dt)
{
  const k1 = a(s).add(new State(s.vx, s.vy, s.vz, 0, 0, 0));
  const k2 = a(s.add(k1.multiply(dt / 2))).add(new State(s.vx + 0.5 * dt * k1.vx, s.vy + 0.5 * dt * k1.vy, s.vz + 0.5 * dt * k1.vz, 0, 0, 0));
  const k3 = a(s.add(k2.multiply(dt / 2))).add(new State(s.vx + 0.5 * dt * k2.vx, s.vy + 0.5 * dt * k2.vy, s.vz + 0.5 * dt * k2.vz, 0, 0, 0));
  const k4 = a(s.add(k3.multiply(dt))).add(new State(s.vx + dt * k3.vx, s.vy + dt * k3.vy, s.vz + dt * k3.vz, 0, 0, 0));
  return s.add(k1.add(k2.multiply(2)).add(k3.multiply(2)).add(k4).multiply(dt / 6));
}

var sunZ;

function transformCoordinates(stateVectors)
{
  // Step 1: Find the geometric center
  let centerX = 0, centerY = 0, centerZ = 0;
  for (let i = 0; i < stateVectors.length; i++)
  {
    centerX += stateVectors[i].x;
    centerY += stateVectors[i].y;
    centerZ += stateVectors[i].z;
  }
  centerX /= stateVectors.length;
  centerY /= stateVectors.length;
  centerZ /= stateVectors.length;

  // Step 2: Designate the new origin O'
  let newOrigin = {x: centerX, y: centerY, z: centerZ};

  // Step 3: Establish the new coordinate system
  let zPrimeAxis = {x: centerX, y: centerY, z: centerZ};
  let zPrimeAxisLength = Math.sqrt(zPrimeAxis.x**2 + zPrimeAxis.y**2 + zPrimeAxis.z**2);
  zPrimeAxis.x /= zPrimeAxisLength;
  zPrimeAxis.y /= zPrimeAxisLength;
  zPrimeAxis.z /= zPrimeAxisLength;

  let xPrimeAxis = {x: -zPrimeAxis.y, y: zPrimeAxis.x, z: 0};
  let xPrimeAxisLength = Math.sqrt(xPrimeAxis.x**2 + xPrimeAxis.y**2 + xPrimeAxis.z**2);
  xPrimeAxis.x /= xPrimeAxisLength;
  xPrimeAxis.y /= xPrimeAxisLength;
  xPrimeAxis.z /= xPrimeAxisLength;

  // Step 4: Express the coordinates of each spacecraft in the new coordinate system
  let transformedCoordinates = [];
  for (let i = 0; i < stateVectors.length; i++) {
    let dx = stateVectors[i].x - newOrigin.x;
    let dy = stateVectors[i].y - newOrigin.y;
    let dz = stateVectors[i].z - newOrigin.z;

    let xPrime = dx * xPrimeAxis.x + dy * xPrimeAxis.y + dz * xPrimeAxis.z;
    let yPrime = dx * (zPrimeAxis.y * xPrimeAxis.z - zPrimeAxis.z * xPrimeAxis.y)
			   + dy * (zPrimeAxis.z * xPrimeAxis.x - zPrimeAxis.x * xPrimeAxis.z)
			   + dz * (zPrimeAxis.x * xPrimeAxis.y - zPrimeAxis.y * xPrimeAxis.x);
    let zPrime = dx * zPrimeAxis.x + dy * zPrimeAxis.y + dz * zPrimeAxis.z;

    transformedCoordinates.push({x: xPrime, y: yPrime, z: zPrime,
			vx: stateVectors[i].vx, vy: stateVectors[i].vy, vz: stateVectors[i].vz,
			r: stateVectors[i].r, g: stateVectors[i].g, b: stateVectors[i].b});
  }

  // Sanity check: Calculate the position of the Sun!
  // Should always be {x:0, y:0, z:-O(AU)}.
  {
    let dx = 0 - newOrigin.x;
    let dy = 0 - newOrigin.y;
    let dz = 0 - newOrigin.z;
    let xPrime = dx * xPrimeAxis.x + dy * xPrimeAxis.y + dz * xPrimeAxis.z;
    let yPrime = dx * (zPrimeAxis.y * xPrimeAxis.z - zPrimeAxis.z * xPrimeAxis.y)
			   + dy * (zPrimeAxis.z * xPrimeAxis.x - zPrimeAxis.x * xPrimeAxis.z)
			   + dz * (zPrimeAxis.x * xPrimeAxis.y - zPrimeAxis.y * xPrimeAxis.x);
    let zPrime = dx * zPrimeAxis.x + dy * zPrimeAxis.y + dz * zPrimeAxis.z;
    if (Math.abs(xPrime) > 1 || Math.abs(yPrime) > 1)
    {
      console.log("The SUN is not where it is supposed to be!");
      console.log({x:xPrime, y:yPrime, z: zPrime});
    }
    sunZ = zPrime;
  }

  return transformedCoordinates;
}

function determinant(r1, r2, r3, r4)
{
	return	r2.x*r3.y*r4.z-r1.x*r3.y*r4.z-r2.y*r3.x*r4.z+r1.y*r3.x*r4.z+
			r1.x*r2.y*r4.z-r1.y*r2.x*r4.z-r2.x*r3.z*r4.y+r1.x*r3.z*r4.y+
			r2.z*r3.x*r4.y-r1.z*r3.x*r4.y-r1.x*r2.z*r4.y+r1.z*r2.x*r4.y+
			r2.y*r3.z*r4.x-r1.y*r3.z*r4.x-r2.z*r3.y*r4.x+r1.z*r3.y*r4.x+
			r1.y*r2.z*r4.x-r1.z*r2.y*r4.x-r1.x*r2.y*r3.z+r1.y*r2.x*r3.z+
			r1.x*r2.z*r3.y-r1.z*r2.x*r3.y-r1.y*r2.z*r3.x+r1.z*r2.y*r3.x;
}

function volume(spacecrafts)
{
	return 1/6*determinant(spacecrafts[0], spacecrafts[1], spacecrafts[2], spacecrafts[3]);
}

function area(s1, s2, s3)
{
    var d = ((s2.z-s3.z)**2+(s2.y-s3.y)**2+(s2.x-s3.x)**2)*((s1.x-s3.x)**2+(s1.z-s3.z)**2+(s1.y-s3.y)**2);
    d -= ((s2.x-s3.x)*(s1.x-s3.x)+(s1.z-s3.z)*(s2.z-s3.z)+(s1.y-s3.y)*(s2.y-s3.y))**2;
    return Math.sqrt(d);
}

function L(r1,r2)
{
  return Math.sqrt((r1.x-r2.x)**2+(r1.y-r2.y)**2+(r1.z-r2.z)**2);
}

function X(u,v)
{
   return {x:u.y*v.z-u.z*v.y, y:u.z*v.x-u.x*v.z, z:u.x*v.y-u.y*v.x};
}

function dot(u,v)
{
  return u.x*v.x+u.y*v.y+u.z*v.z;
}

function norm(R)
{
	return Math.sqrt(R.x*R.x + R.y*R.y + R.z*R.z);
}

function smul(s, v)
{
  return {x:s*v.x, y:s*v.y, z:s*v.z};
}

function Mmul(M, v)
{
  return {x:M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z,
          y:M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z,
          z:M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z};
}

function MP(M, K)
{
  let R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
  {
    R[i][j] = M[i][j] + K[i][j];
  }
  return R;
}

function MM(M, K)
{
  let R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
  {
    R[i][j] += M[i][k]*K[k][j];
  }
  return R;
}

function sM(s, M)
{
  let R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
  {
    R[i][j] = s * M[i][j];
  }
  return R;
}

function vadd(u, v)
{
  return {x:u.x+v.x, y:u.y+v.y, z:u.z+v.z};
}

function inv3x3M(M)
{
  const a = M[0][0], b = M[0][1], c = M[0][2];
  const d = M[1][0], e = M[1][1], f = M[1][2];
  const g = M[2][0], h = M[2][1], i = M[2][2];

  const D = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);

  return [[(e*i - f*h) / D, (c*h - b*i) / D, (b*f - c*e) / D],
          [(f*g - d*i) / D, (a*i - c*g) / D, (c*d - a*f) / D],
          [(d*h - e*g) / D, (b*g - a*h) / D, (a*e - b*d) / D]];
}

class TS
{
  constructor(size)
  {
    this.size = size;
    this.data = [];
  }
  push(element)
  {
    this.data.push(element);
    if (this.data.length > this.size) this.data.shift();
  }
}

rtsM = new TS(3);
rtsT = new TS(3);
rtsE = new TS(3);
rtsW = new TS(3);

function trT(states, DT)
{
  let r12 = {x:states[1].x - states[0].x, y:states[1].y - states[0].y, z: states[1].z - states[0].z};
  let r13 = {x:states[2].x - states[0].x, y:states[2].y - states[0].y, z: states[2].z - states[0].z};
  let r14 = {x:states[3].x - states[0].x, y:states[3].y - states[0].y, z: states[3].z - states[0].z};

  rtsT.push({r12:r12, r13:r13, r14:r14});

  return doTR(rtsT, DT);
}

function trM(states, DT)
{
  let a = L(states[0], states[1]);
  let b = L(states[0], states[2]);
  let c = L(states[1], states[2]);
  let d = L(states[0], states[3]);
  let e = L(states[1], states[3]);
  let f = L(states[2], states[3]);

  let r12={x:a, y:0, z:0};
  let r13={x:(a**2 + b**2 - c**2)/(2*a), y:Math.sqrt(b**2 - ((a**2 + b**2 - c**2)/(2*a))**2), z:0};
  let r14={x:(a**2 + d**2 - e**2)/(2*a), y:(b**2 + d**2 - f**2 - ((a**2 + d**2 - e**2)*(a**2 + b**2 - c**2))/(2*a**2))/(2*Math.sqrt(b**2 - ((a**2 + b**2 - c**2)/(2*a))**2)), z:Math.sqrt(d**2 - ((a**2 + d**2 - e**2)/(2*a))**2 - ((b**2 + d**2 - f**2 - ((a**2 + d**2 - e**2)*(a**2 + b**2 - c**2))/(2*a**2))/(2*Math.sqrt(b**2 - ((a**2 + b**2 - c**2)/(2*a))**2)))**2)};

  rtsM.push({r12:r12, r13:r13, r14:r14});

  return doTR(rtsM, DT);
}

function trE(states, DT)
{
  // What we want is the same calculation with respect to the "satellite" frame.

  let r12 = {x:states[1].x - states[0].x, y:states[1].y - states[0].y, z: states[1].z - states[0].z};
  let r13 = {x:states[2].x - states[0].x, y:states[2].y - states[0].y, z: states[2].z - states[0].z};
  let r14 = {x:states[3].x - states[0].x, y:states[3].y - states[0].y, z: states[3].z - states[0].z};

  let ex = states.reduce((a,b)=>({x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}));
  let nx = norm(ex);
  ex.x /= nx;
  ex.y /= nx;
  ex.z /= nx;

  let ev = states.reduce((a,b)=>({vx:a.vx+b.vx,vy:a.vy+b.vy,vz:a.vz+b.vz}));

  let ez = {x:ex.y*ev.vz - ex.z*ev.vy, y:ex.z*ev.vx - ex.x*ev.vz, z:ex.x*ev.vy - ex.y*ev.vx};
  let nz = norm(ez);
  ez.x /= nz;
  ez.y /= nz;
  ez.z /= nz;

  let ey = {x:ez.y*ex.z-ez.z*ex.y, y:ez.z*ex.x -ez.x*ex.z, z:ez.x*ex.y - ez.y*ex.x};

  let R12 = {x:dot(r12,ex), y:dot(r12,ey), z:dot(r12,ez)};
  let R13 = {x:dot(r13,ex), y:dot(r13,ey), z:dot(r13,ez)};
  let R14 = {x:dot(r14,ex), y:dot(r14,ey), z:dot(r14,ez)};

  rtsE.push({r12:R12, r13:R13, r14:R14});

  return doTR(rtsE, DT);
}

function Sagnac(ri, rj, vi, vj)
{
  const c = 299.792458; // Mm/s

  let t1i = (dot(vi,ri) + Math.sqrt(dot(vi,ri)**2 + (c**2-dot(vi,vi))*dot(ri,ri))) / (c**2 - dot(vi,vi));
  let rij = vadd(rj, smul(-1, ri));
  let vij = vadd(vj, smul(-1, vi));
  let bij = dot(vadd(rij, smul(t1i, vij)), vj)
  let tij = (bij + Math.sqrt(bij**2 + (c**2 - dot(vj,vj))*dot(vadd(rij, smul(t1i, vij)),vadd(rij, smul(t1i, vij))))) / (c**2 - dot(vj,vj));
  let tj1 = Math.sqrt(dot(vadd(rj, smul((t1i+tij), vj)), vadd(rj, smul((t1i+tij), vj)))) / c;

  return t1i + tij + tj1;
}

function trW(DT)
{
  if (rtsT.data.length != 3) return NaN;

  // Midpoint position
  let r12 = rtsT.data[1].r12;
  let r13 = rtsT.data[1].r13;
  let r14 = rtsT.data[1].r14;

  // Velocities of 234 wrt. satellite 1
  let v12 = smul(1/(2*DT), vadd(rtsT.data[2].r12, smul(-1, rtsT.data[0].r12)));
  let v13 = smul(1/(2*DT), vadd(rtsT.data[2].r13, smul(-1, rtsT.data[0].r13)));
  let v14 = smul(1/(2*DT), vadd(rtsT.data[2].r14, smul(-1, rtsT.data[0].r14)));

  let w123 = Sagnac(r12, r13, v12, v13) - Sagnac(r13, r12, v13, v12);
  let w134 = Sagnac(r13, r14, v13, v14) - Sagnac(r14, r13, v14, v13);
  let w142 = Sagnac(r14, r12, v14, v12) - Sagnac(r12, r14, v12, v14);

//console.log("w134: " + w134 + ", w142: " + w142 + ", w123: " + w123);

  //
  //
  // From this point on, all we are allowed to use are w1..w3 and rtsM
  //
  //

  r12 = rtsM.data[1].r12;
  r13 = rtsM.data[1].r13;
  r14 = rtsM.data[1].r14;

  // Velocities of 234 wrt. satellite 1
  v12 = smul(1/(2*DT), vadd(rtsM.data[2].r12, smul(-1, rtsM.data[0].r12)));
  v13 = smul(1/(2*DT), vadd(rtsM.data[2].r13, smul(-1, rtsM.data[0].r13)));
  v14 = smul(1/(2*DT), vadd(rtsM.data[2].r14, smul(-1, rtsM.data[0].r14)));

  let w = {x:0, y:0, z:0};
  let eps = 1e-10;
  let dx = {x:eps, y:0, z:0};
  let dy = {x:0, y:eps, z:0};
  let dz = {x:0, y:0, z:eps};

  function m1234(w, r12, r13, r14, v12, v13, v14)
  {
    let w12 = vadd(v12, X(r12, w));
    let w13 = vadd(v13, X(r13, w));
    let w14 = vadd(v14, X(r14, w));

    let m =
    {
      m123: Sagnac(r12, r13, w12, w13) - Sagnac(r13, r12, w13, w12),
      m134: Sagnac(r13, r14, w13, w14) - Sagnac(r14, r13, w14, w13),
      m142: Sagnac(r14, r12, w14, w12) - Sagnac(r12, r14, w12, w14)
    };
    return m;
  };

  let count = 50;
  let m;
  let D = {x:0, y:0, z:0};
  let F = {x:0, y:0, z:0};
  do
  {
    w = vadd(w, D);

    m = m1234(w, r12, r13, r14, v12, v13, v14);

    let mx = m1234(vadd(w, dx), r12, r13, r14, v12, v13, v14);
    let my = m1234(vadd(w, dy), r12, r13, r14, v12, v13, v14);
    let mz = m1234(vadd(w, dz), r12, r13, r14, v12, v13, v14);

    let J = [[(mx.m123-m.m123)/eps, (my.m123-m.m123)/eps, (mz.m123-m.m123)/eps],
             [(mx.m134-m.m134)/eps, (my.m134-m.m134)/eps, (mz.m134-m.m134)/eps],
             [(mx.m142-m.m142)/eps, (my.m142-m.m142)/eps, (mz.m142-m.m142)/eps]];

    F = {x:w123 - m.m123, y:w134 - m.m134, z:w142 - m.m142};
    D = Mmul(inv3x3M(J), F);
  
  } while (--count > 0 && Math.sqrt(dot(F,F)) > 1e-17);

//console.log("COUNT = " + count);
//console.log(D);
//console.log("m134: " + m.m134 + ", m142: " + m.m142 + ", m123: " + m.m123);
//console.log(w);

  // We can now calculate the velocities of the 234 vertices wrt. 1
//  let u12 = X(w, smul(0.5, vadd(rtsM.data[0].r12, rtsM.data[1].r12)));
//  let u13 = X(w, smul(0.5, vadd(rtsM.data[0].r13, rtsM.data[1].r13)));
//  let u14 = X(w, smul(0.5, vadd(rtsM.data[0].r14, rtsM.data[1].r14)));
//console.log(u12);
//console.log(u13);
//console.log(u14);

  let W = Math.sqrt(dot(w,w));
  let theta = W*DT;
  let I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  let K = [[0, -w.z/W, w.y/W], [w.z/W, 0, -w.x/W], [-w.y/W, w.x/W, 0]];
  let R = MP(I, MP(sM(Math.sin(theta), K), sM(1-Math.cos(theta), MM(K,K))));

  // and adjust their relative positions, removing the fictitious force
  r12 = Mmul(R, rtsM.data[0].r12);
  r13 = Mmul(R, rtsM.data[0].r13);
  r14 = Mmul(R, rtsM.data[0].r14);

  rtsW.push({r12:r12, r13:r13, r14:r14});

  // The middle position is left alone
  rtsW.push(rtsM.data[1]);

  theta = -W*DT;
  R = MP(I, MP(sM(Math.sin(theta), K), sM(1-Math.cos(theta), MM(K,K))));

  r12 = Mmul(R, rtsM.data[2].r12);
  r13 = Mmul(R, rtsM.data[2].r13);
  r14 = Mmul(R, rtsM.data[2].r14);

  rtsW.push({r12:r12, r13:r13, r14:r14});

  return doTR(rtsW, DT);
}

function doTR(rts, DT)
{
  if (rts.data.length == 3)
  {
    let a12={x:(rts.data[0].r12.x+rts.data[2].r12.x-2*rts.data[1].r12.x)/DT**2,
             y:(rts.data[0].r12.y+rts.data[2].r12.y-2*rts.data[1].r12.y)/DT**2,
             z:(rts.data[0].r12.z+rts.data[2].r12.z-2*rts.data[1].r12.z)/DT**2};
    let a13={x:(rts.data[0].r13.x+rts.data[2].r13.x-2*rts.data[1].r13.x)/DT**2,
             y:(rts.data[0].r13.y+rts.data[2].r13.y-2*rts.data[1].r13.y)/DT**2,
             z:(rts.data[0].r13.z+rts.data[2].r13.z-2*rts.data[1].r13.z)/DT**2};
    let a14={x:(rts.data[0].r14.x+rts.data[2].r14.x-2*rts.data[1].r14.x)/DT**2,
             y:(rts.data[0].r14.y+rts.data[2].r14.y-2*rts.data[1].r14.y)/DT**2,
             z:(rts.data[0].r14.z+rts.data[2].r14.z-2*rts.data[1].r14.z)/DT**2};
    let r12=rts.data[1].r12;
    let r13=rts.data[1].r13;
    let r14=rts.data[1].r14;

/*  // Alternative method: explicitly calculating the trace
    let A = [[a12.x, a13.x, a14.x], [a12.y, a13.y, a14.y], [a12.z, a13.z, a14.z]];
    let R = [[r12.x, r13.x, r14.x], [r12.y, r13.y, r14.y], [r12.z, r13.z, r14.z]];
    let T = MM(A, inv3x3M(R));
    console.log(T);
    let tr = T[0][0] + T[1][1] + T[2][2];
*/

    let V2 = dot(r12,X(r13,r14));
    let V3 = dot(r13,X(r14,r12));
    let V4 = dot(r14,X(r12,r13));

    let tr = dot(a12,X(r13,r14))/V2 + dot(a13,X(r14,r12))/V3 + dot(a14,X(r12,r13))/V4;

    return tr;
  }
  return NaN;
}

function eccentricity(state)
{
  let r = Math.sqrt(state.x*state.x + state.y*state.y + state.z*state.z);
  let v = Math.sqrt(state.vx*state.vx + state.vy*state.vy + state.vz*state.vz);
  let h = Math.sqrt(GM*r);
  let e = (v*v/GM - 1/r)/(v*v/GM + 1/r);
  return e; 
}

function avgeccentricity(states)
{
  let x = 0; let y = 0; let z = 0;
  let vx = 0; let vy = 0; let vz = 0;
  for (let i=0; i<states.length; i++)
  {
    x += states[i].x;
    y += states[i].y;
    z += states[i].z;
    vx += states[i].vx;
    vy += states[i].vy;
    vz += states[i].vz;
  }
  x /= states.length; y /= states.length; z /= states.length;
  vx /= states.length; vy /= states.length; vz /= states.length;

  let hx = y*vz - z*vy;
  let hy = z*vx - x*vz;
  let hz = x*vy - y*vx;
  let h = Math.sqrt(hx*hx + hy*hy + hz*hz);
  let r = Math.sqrt(x*x + y*y + z*z);
  let v = Math.sqrt(vx*vx + vy*vy + vz*vz);
//  let a = 1 / (2/r - v*v/GM);
//  let e = Math.sqrt(1 - h*h/a/GM);

  let e = Math.sqrt(1 + h*h/GM/GM*(v*v-2*GM/r));
  return e; 
}

function orbitalPeriod(stateVector)
{
  let r = Math.sqrt(stateVector.x*stateVector.x + stateVector.y*stateVector.y + stateVector.z*stateVector.z);
  let v = Math.sqrt(stateVector.vx*stateVector.vx + stateVector.vy*stateVector.vy + stateVector.vz*stateVector.vz);
  let epsilon = 0.5*v*v - GM/r;
  let a = -GM/2/epsilon;
  let T = 2*Math.PI*Math.sqrt(a*a*a/GM);
  return T; 
}

function meanOrbitalPeriod(stateVectors)
{
  let sumX = 0; 
  let sumY = 0;
  let sumZ = 0;
  let sumVx = 0;
  let sumVy = 0;
  let sumVz = 0;
  for (let i=0; i<stateVectors.length; i++) {
    sumX += stateVectors[i].x;
    sumY += stateVectors[i].y;
    sumZ += stateVectors[i].z;
    sumVx += stateVectors[i].vx;
    sumVy += stateVectors[i].vy;
    sumVz += stateVectors[i].vz;
  }
  let mean = {
    x: sumX/stateVectors.length,
    y: sumY/stateVectors.length, 
    z: sumZ/stateVectors.length,
    vx: sumVx/stateVectors.length,
    vy: sumVy/stateVectors.length, 
    vz: sumVz/stateVectors.length
  };
  return orbitalPeriod(mean); 
}

var time = 0;
var inPROP = false;

function propagate(dontMove = false, forward = true)
{
	if (inPROP) return;
	inPROP = true;

	const DT = forward ? 1800 : -1800;
	const NT = 48;
	var traceM = 0;
	var traceT = 0;
	var traceE = 0;
	var traceW = 0;

	if (!dontMove)
	{
	for (let i = 0; i < NT; i++)
	{
		for (let j = 0; j < states.length; j++)
		{
			states[j] = rk4(states[j], DT);
		}
		traceM = trM(states, DT);
		traceT = trT(states, DT);
		traceE = trE(states, DT);
        traceW = trW(DT);
		time += DT / 86400.0;
	}
	}

	let theVolume = volume(transformCoordinates(states)).toFixed(2);
	if (theVolume[0] != '-') theVolume = "&nbsp;" + theVolume;
	traceM = traceM.toPrecision(5);
	traceT = traceT.toPrecision(5);
	traceE = traceE.toPrecision(5);
	traceW = traceW.toPrecision(5);
	if (traceM[0] != '-') traceM = "&nbsp;" + traceM;
	if (traceT[0] != '-') traceT = "&nbsp;" + traceT;
	if (traceE[0] != '-') traceE = "&nbsp;" + traceE;
	if (traceW[0] != '-') traceW = "&nbsp;" + traceW;

	document.getElementById("volume").innerHTML = theVolume;
	document.getElementById("trM").innerHTML = traceM;
	document.getElementById("trT").innerHTML = traceT;
	document.getElementById("trE").innerHTML = traceE;
	document.getElementById("trW").innerHTML = traceW;
	document.getElementById("time").innerText = time.toFixed(2);

    //window.requestAnimationFrame(render);
	render();

	inPROP = false;
}

function drawShinySphere(ctx, x, y, radius, r, g, b)
{
  // Create a radial gradient
  const gradient = ctx.createRadialGradient(x, y, 0, x, y, radius);

  // Add color stops to create the shiny effect
  gradient.addColorStop(0, `rgba(${r}, ${g}, ${b}, 1)`); // Bright center
  gradient.addColorStop(0.8, `rgba(${Math.round(r/2.5)}, ${Math.round(g/2.5)}, ${Math.round(b/2.5)}, 0.8)`); // Slightly darker
  gradient.addColorStop(1, 'rgba(0, 0, 0, 0)'); // Faded edges

  // Set the fill style to the radial gradient
  ctx.fillStyle = gradient;

  // Draw the circle
  ctx.beginPath();
  ctx.arc(x, y, radius, 0, 2 * Math.PI);
  ctx.closePath();
  ctx.fill();
}

function setView(v)
{
	var v = document.getElementById("view").value;
	var d = 3 + 0.27 * document.getElementById("zoom").value;

	if (v == 1) { camera.x = d; camera.y = 0; camera.z = 0; }
	if (v == 2) { camera.x = 0; camera.y = d; camera.z = 0; }
	if (v == 3) { camera.x = 0; camera.y = 0; camera.z = d; }
	camera.phi = camera.theta = 0;
	view = v;
	if (timer == 0) render();
}

function render()
{
	var spacecrafts = transformCoordinates(states);
	var v = document.getElementById("view").value * 1;

	spacecrafts.forEach(sc =>
	{
		var tempX, tempY, tempZ;

		switch (v)
		{
		case 1:
			tempX = sc.x * Math.cos(camera.phi) + sc.y * Math.sin(camera.phi);
			sc.y = -sc.x * Math.sin(camera.phi) + sc.y * Math.cos(camera.phi);
			sc.z = sc.z * Math.cos(camera.theta) + tempX * Math.sin(camera.theta);
			sc.x = tempX;
			break;
		case 2:
			tempY = sc.y * Math.cos(camera.theta) - sc.x * Math.sin(camera.theta);
			sc.x = sc.y * Math.sin(camera.theta) + sc.x * Math.cos(camera.theta);
			sc.z = sc.z * Math.cos(camera.phi) - tempY * Math.sin(camera.phi);
			sc.y = tempY;
			break;
		case 3:
			tempZ = sc.z * Math.cos(camera.theta) - sc.y * Math.sin(camera.theta);
			sc.y = sc.z * Math.sin(camera.theta) + sc.y * Math.cos(camera.theta);
			sc.x = sc.x * Math.cos(camera.phi) - tempZ * Math.sin(camera.phi);
			sc.z = tempZ;
			break;
		}
	});

    // Clear the canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    if (view == 3)
    {
	  // Draw sun
	  ctx.fillStyle = 'yellow';
	  ctx.beginPath();
      let cw = Math.min(canvas.width, canvas.height);
	  ctx.arc(canvas.width/2 + camera.phi*4/Math.PI*cw/2, canvas.height/2 + camera.theta*4/Math.PI*cw/2, -5*AU/sunZ, 0, 2 * Math.PI);
	  ctx.fill();
    }

	// Sort the indices array based on the 'z' property of the corresponding objects in the original array

	if (view == 1) indices.sort((a, b) => spacecrafts[a].x - spacecrafts[b].x);
	if (view == 2) indices.sort((a, b) => spacecrafts[a].y - spacecrafts[b].y);
	if (view == 3) indices.sort((a, b) => spacecrafts[a].z - spacecrafts[b].z);

	var views = [];

	// Iterate through the sorted indices array and access the original array elements
	indices.forEach(index =>
	{
	  var spacecraft = spacecrafts[index];

      // Calculate distance from the camera
      const dx = spacecraft.x - camera.x;
      const dy = spacecraft.y - camera.y;
      const dz = spacecraft.z - camera.z;

      const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

      // Calculate sphere size based on distance
      const size = 500 / Math.pow(distance, 1.5);

      // Calculate screen coordinates
      var screenX, screenY;
	  var tempX, tempY, tempZ;

	  const canvSize = Math.min(canvas.width, canvas.height);

      if (view == 1)
      {
		if (dx > 0) return;
        tempX = (spacecraft.y / (camera.x - spacecraft.x)) * canvSize;
        tempY = (spacecraft.z / (camera.x - spacecraft.x)) * canvSize;
        //tempZ = (spacecraft.x / (camera.x - spacecraft.x)) * canvSize;
      }
      if (view == 2)
      {
		if (dy > 0) return;
        tempX = (spacecraft.z / (camera.y - spacecraft.y)) * canvSize;
        tempY = (spacecraft.x / (camera.y - spacecraft.y)) * canvSize;
        //tempZ = (spacecraft.y / (camera.y - spacecraft.y)) * canvSize;
      }
      if (view == 3)
      {
		if (dz > 0) return;
        tempX = (spacecraft.x / (camera.z - spacecraft.z)) * canvSize;
        tempY = (spacecraft.y / (camera.z - spacecraft.z)) * canvSize;
        //tempZ = (spacecraft.z / (camera.z - spacecraft.z)) * canvSize;
      }

	  //screenX = tempX * Math.cos(camera.phi);
	  //screenY = tempY * Math.cos(camera.theta) - Math.sqrt(tempX**2 + tempY**2) * Math.sin(camera.theta);

	  //screenX = canvas.width / 2 + screenX;
	  //screenY = canvas.height / 2 - screenY;
	  screenX = canvas.width / 2 + tempX;
	  screenY = canvas.height / 2 - tempY;


      // Draw sphere
	  // drawShinySphere(ctx, screenX, screenY, size, spacecraft.r, spacecraft.g, spacecraft.b);
	  views.push({X:screenX, Y:screenY, S:size, R:spacecraft.r, G:spacecraft.g, B:spacecraft.b});
    });

  // Draw lines connecting pairs of spacecraft
  for (let i = 0; i < views.length; i++)
  {
    const spacecraft1 = views[i];
    drawShinySphere(ctx, spacecraft1.X, spacecraft1.Y, spacecraft1.S, spacecraft1.R, spacecraft1.G, spacecraft1.B);
    for (let j = i + 1; j < views.length; j++)
    {
      const spacecraft2 = views[j];

      // Create a gradient for the line
      const gradient = ctx.createLinearGradient(spacecraft1.X, spacecraft1.Y, spacecraft2.X, spacecraft2.Y);
      gradient.addColorStop(0, `rgba(${spacecraft1.R}, ${spacecraft1.G}, ${spacecraft1.B}, ${Math.pow(spacecraft1.S,3) / 5000})`);
      gradient.addColorStop(1, `rgba(${spacecraft2.R}, ${spacecraft2.G}, ${spacecraft2.B}, ${Math.pow(spacecraft2.S,3) / 5000})`);

      // Draw the line
      ctx.beginPath();
      ctx.moveTo(spacecraft1.X, spacecraft1.Y);
      ctx.lineTo(spacecraft2.X, spacecraft2.Y);
      ctx.strokeStyle = gradient;
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  topView.update(states);

  document.getElementById('phi').innerText = (-camera.phi*180/Math.PI).toFixed(2);
  document.getElementById('theta').innerText = (camera.theta*180/Math.PI).toFixed(2);
}

function zoom()
{
	var d = 3 + 0.27 * document.getElementById("zoom").value;
	if (view == 1) camera.x = d;
	if (view == 2) camera.y = d;
	if (view == 3) camera.z = d;
	if (timer == 0) render();
}

var timer = 0;

function start()
{
	document.getElementById('skipForward').disabled = true;
	document.getElementById('skipBack').disabled = true;
	if (timer == 0) timer = setInterval(propagate, 50);
}

function stop()
{
	document.getElementById('skipForward').disabled = false;
	document.getElementById('skipBack').disabled = false;
	if (timer != 0)
	{
		clearInterval(timer);
		timer = 0;
	}
}

function play()
{
	document.getElementById("init").disabled = true;
	if (timer == 0) start();
	else stop();
}

function skipForward()
{
	document.getElementById("init").disabled = true;
  if (timer == 0) propagate(false, true);
}

function skipBack()
{
	document.getElementById("init").disabled = true;
  if (timer == 0) propagate(false, false);
}

function onProp()
{
	if (timer == 0) propagate(true);
}

var canvas;
var ctx;

function onInit(init)
{
	if (init >= 0 && init < spacecraftsets.length)
	{
		spacecrafts = spacecraftsets[init].states;
		topView.bufferSize = spacecraftsets[init].tail;
	}

	states = spacecrafts.map(spacecraft => new State(spacecraft.x, spacecraft.y, spacecraft.z,
													spacecraft.vx, spacecraft.vy, spacecraft.vz,
													spacecraft.r, spacecraft.g, spacecraft.b));

	// Create an array of indices
	indices = spacecrafts.map((_, index) => index);

	document.getElementById("e").innerText = avgeccentricity(states).toFixed(4);
	document.getElementById("T").innerText = (meanOrbitalPeriod(states)/86400).toFixed(2);
	propagate(true);
}

window.addEventListener('DOMContentLoaded', () =>
{
	const params = new URLSearchParams(window.location.search);
	if (params.get('view')) view = params.get('view');
	if (params.get('init')) init = params.get('init');
	if (params.get('w')) MOD.w = params.get('w');
	if (params.get('m')) MOD.m = params.get('m');
	if (params.get('y')) MOD.y = params.get('y');

	const select = document.getElementById("init");
	for (let i = 0; i < spacecraftsets.length; i++)
	{
		const option = document.createElement('option');
		option.value = i;
		option.text = spacecraftsets[i].name;
		select.appendChild(option);
	}
	select.value = init;

	doResize();
	onInit(init);
	document.getElementById("view").value = view;
	setView(view);
});


function doResize()
{
	const body = document.querySelector('body');
    var ratioX = window.innerWidth / 1200.0;
    var ratioY = window.innerHeight / 800.0;
	var ratio = Math.min(ratioX, ratioY);
	if (ratio < 0.4) ratio = 0.4;
	if (ratio > 3) ratio = 3;
    body.style.zoom = ratio;

	const visualization = document.getElementById('visualization');
	visualization.innerHTML = "";
	canvas = document.createElement('canvas');
	ctx = canvas.getContext('2d');
	canvas.width = visualization.clientWidth;
	canvas.height = visualization.clientHeight;
	visualization.appendChild(canvas);

	let z = ([camera.x, camera.y, camera.z])[view - 1];
	document.getElementById("zoom").value=(z - 3)*100/27;

	function doZoom(delta, step=5)
	{
		var s = document.getElementById("zoom");
		var z = s.value;
		z = z - z%step;
		if (delta < 0 && z > 0) z -= step;
		if (delta > 0 && z < 100) z += step;
		s.value = z;
		zoom();
	}

	// Add a scroll event listener to the canvas
	canvas.addEventListener("wheel", function(event)
	{
		event.preventDefault();

		doZoom(event.deltaY, 5);
	});

	function doClick(e)
	{
		if (!canvas.wasdragging && !canvas.wasscaling) play();
		canvas.wasdragging = false;
		canvas.wasscaling = false;
	}

	function distance(touches)
	{
		return Math.sqrt((touches[1].clientX - touches[0].clientX)**2
						+ (touches[1].clientY - touches[0].clientY)**2);
	}

	canvas.addEventListener("click", doClick);

	function doStart(e)
	{
		canvas.dragging = true;
		if (e.touches)
		{
			canvas.lastX = e.touches[0].clientX;
			canvas.lastY = e.touches[0].clientY;
			if (e.touches.length === 2)
			{
				canvas.lastD = distance(e.touches);
			}
		}
		else
		{
		canvas.lastX = e.offsetX;
		canvas.lastY = e.offsetY;
		}
	}

	canvas.addEventListener("mousedown", doStart);
	canvas.addEventListener("touchstart", (e) =>
	{
		e.preventDefault();
		doStart(e);
	});

	function doMove(e)
	{
		if (canvas.dragging)
		{
			let dx = 0;
			let dy = 0;

			if (e.touches)
			{
				if (e.touches.length === 2)
				{
					let d = distance(e.touches);
					if (canvas.lastD != d)
					{
						doZoom(canvas.lastD - d, 1);
						canvas.wasscaling = true;
						canvas.lastD = d;
					}
				}
				else
				{
					dx = e.touches[0].clientX - canvas.lastX;
					dy = e.touches[0].clientY - canvas.lastY;

					canvas.lastX = e.touches[0].clientX;
					canvas.lastY = e.touches[0].clientY;
					if (dx != 0 || dy != 0) canvas.wasdragging = true;
				}
			}
			else
			{
				dx = e.offsetX - canvas.lastX;
				dy = e.offsetY - canvas.lastY;

				canvas.lastX = e.offsetX;
				canvas.lastY = e.offsetY;
				canvas.wasdragging = true;
			}

			camera.phi -= .5 * dx / 180 * Math.PI;
			camera.theta -= .5 * dy / 180 * Math.PI;

			if (camera.phi < -0.5*Math.PI) camera.phi = -0.5*Math.PI;
			if (camera.phi > 0.5*Math.PI) camera.phi = 0.5*Math.PI;
			if (camera.theta < -0.5*Math.PI) camera.theta = -0.5*Math.PI;
			if (camera.theta > 0.5*Math.PI) camera.theta = 0.5*Math.PI;
			if (!timer) render();
		}
	}

	canvas.addEventListener("mousemove", doMove);
	canvas.addEventListener("touchmove", (e) =>
	{
		e.preventDefault();
		doMove(e);
	});

	function doStop(e)
	{
		canvas.dragging = false;
	}

	canvas.addEventListener("mouseup", doStop);
	canvas.addEventListener("touchend", (e) =>
	{
		//const touchcancel = new TouchEvent('touchcancel',
		//{
		//	bubbles: true,
		//	cancelable: true
		//});
		//e.target.dispatchEvent(touchcancel);
		doStop();
		doClick();
	});

	topView.init(canvas, canvas.width-topView.width, canvas.height-topView.height);
}

let topView =
{
  width: 200,
  height: 200,
  context: null,
  buffer: [],
  bufferSize: 300,
  x: 0, // x coordinate of top left of view
  y: 0, // y coordinate of top left of view

  init(canvas, x, y)
  {
    this.context = canvas.getContext('2d');
    this.x = x;
    this.y = y;

  },

  update(stateVectors)
  {
    // Average state vectors
    let avgX = 0;
    let avgY = 0;
    let avgZ = 0;
    stateVectors.forEach(v =>
    {
      avgX += v.x;
      avgY += v.y;
      avgZ += v.z;
    });
    avgX /= stateVectors.length;
    avgY /= stateVectors.length;
    avgZ /= stateVectors.length;
    let avgState = {x: avgX / AU, y: -avgY / AU, z: avgZ / AU};

    // Add latest state to buffer
	if (timer != 0)
	{
    	this.buffer.push(avgState);
    	if (this.buffer.length > this.bufferSize)
    	{
    	  this.buffer.shift();
    	}
	}

    // Clear view area only
    this.context.clearRect(this.x, this.y, this.width, this.height);

	// Draw sun
	this.context.fillStyle = 'yellow';
	this.context.beginPath();
	this.context.arc(this.x+this.width/2, this.y+this.height/2, 8, 0, 2 * Math.PI);
	this.context.fill();

    // Scale state vector coordinates to view
    let x = (avgState.x + 2.5) * this.width / 5 + this.x;
    let y = (avgState.y + 2.5) * this.height / 5 + this.y;

    // Draw current position
    this.context.fillStyle = '#8FF';
	this.context.beginPath();
	this.context.arc(x, y, 4, 0, 2 * Math.PI);
	this.context.fill();

	// Draw trail
	this.context.lineWidth = 1;
	this.context.beginPath();
	this.context.moveTo(x, y);
	for (let i = this.buffer.length - 2; i >= 0; i -= 5)
	{
		let state = this.buffer[i];
		let x = (state.x + 2.5) * this.width / 5 + this.x;
		let y = (state.y + 2.5) * this.height / 5 + this.y;
		this.context.lineTo(x, y);
		// Fade color
		this.context.strokeStyle = `rgba(224, 224, 255, ${(i/this.buffer.length)**0.125})`;
		this.context.stroke();
		this.context.beginPath();
		this.context.moveTo(x, y);
	}//);

    // Draw border around view area
	//this.context.beginPath();
    this.context.strokeStyle = 'lightgray';
    this.context.strokeRect(this.x, this.y, this.width, this.height);

    this.context.stroke();
  }
}

window.addEventListener('resize', () => { doResize(); onProp(); });
