const int MAX_MARCHING_STEPS = 127;
const float MIN_DIST = 0.0;
const float MAX_DIST = 10.0;
const float EPSILON = 0.0001;
const vec4 ORIGIN = vec4(0,0,0,1);

const float halfIdealCubeWidthKlein = 0.5773502692;
const vec4 idealCubeCornerKlein = vec4(halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, 1.0);

/// for {4,3,6}
 const float sphereRad = 1.0;
 const float horosphereSize = 2.6;  //spheres intersect

// const float sphereRad = 0.5;
// const float horosphereSize = 1.5;  // horospheres intersect

/// for {4,3,7}
//const float sphereRad = 1.0;      //spheres intersect
//const float planeOffset = 0.75;

uniform vec2 screenResolution;
uniform vec3 cameraPos;
uniform vec4 cameraQuat;
uniform float fov;
uniform mat4 generators[6];
uniform mat4 invGenerators[6];
uniform mat4 currentBoost;
uniform int maxSteps;
uniform float halfCubeWidthKlein;

//-------------------------------------------------------
//Hyperboloid Functions
//-------------------------------------------------------

float lorentzDot(vec4 u, vec4 v){
  return u.w*v.w - u.x*v.x - u.y*v.y - u.z*v.z;
} /// on hyperbolold if lorentzDot(u,u) = 1, so w*w = 1 + x*x + y*y + z*z

vec4 projectToHyperboloid(vec4 v) {  //projects to either hyperboloid depending if input
  //is timelike or spacelike
  return v/sqrt(abs(lorentzDot(v,v)));
}

vec4 projectToKlein(vec4 v){
  return v/v.w;
}

vec3 qtransform( vec4 q, vec3 v ){
  return v + 2.0*cross(cross(v, -q.xyz ) + q.w*v, -q.xyz);
}

vec4 qinverse(vec4 q){
  float n = q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w;
  return vec4(-q.x/n, -q.y/n, -q.z/n, q.w/n);
}

//-------------------------------------------------------
//Hyperbolic Math functions
//-------------------------------------------------------
float cosh(float x){
  float eX = exp(x);
  return (0.5 * (eX + 1.0/eX));
}
float sinh(float x){
  float eX = exp(x);
  return (0.5 * (eX - 1.0/eX));
}
float acosh(float x){ //must be more than 1
  return log(x + sqrt(x*x-1.0));
}
float asinh(float x){
  return log(x + sqrt(x*x+1.0));
}

float hypNorm(vec4 v){
  return sqrt(abs(lorentzDot(v,v)));
}

vec4 lorentzNormalize(vec4 v){  // cannot do to a light like vector
  return v/hypNorm(v);  // note that this is the same function as projectToHyperboloid
}

float hypDistance(vec4 u, vec4 v){
  float bUV = lorentzDot(u,v);
  return acosh(bUV);
}

vec4 vPrimeFromV(vec4 u, vec4 v){  // given points u and v on hyperboloid, make
  // the point vPrime for use in parametrising the geodesic from u through v
  vec4 w = v - lorentzDot(u, v)*u;
  return (1.0/hypNorm(w)*w);
}

vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist){ // get point on
  // hyperboloid at distance dist on the geodesic from u through v
  return u*cosh(dist) + vPrime*sinh(dist);
}

vec4 tangentVectorOnGeodesic(vec4 u, vec4 vPrime, float dist){
  // note that this point has lorentzDot with itself of -1, so it is on other hyperboloid
  return u*sinh(dist) + vPrime*cosh(dist);
}

vec4 pointOnGeodesicAtInfinity(vec4 u, vec4 vPrime){ // returns point on the light
  // cone intersect Klein model corresponding to the point at infinity on the
  // geodesic through u and v
  return projectToKlein(u + vPrime);
}

mat4 translateByVector(vec3 v) { // trickery from Jeff Weeks' Curved Spaces app
  float dx = v.x;
  float dy = v.y;
  float dz = v.z;
  float len = sqrt(dx*dx + dy*dy + dz*dz);
  if (len == 0.0){
    return mat4(1.0);
  }
  else{
      dx /= len;
      dy /= len;
      dz /= len;
      mat4 m = mat4(vec4(0, 0, 0, dx),
                    vec4(0, 0, 0, dy),
                    vec4(0, 0, 0, dz),
                    vec4(dx,dy,dz, 0));
      mat4 m2 = m*m;
      float c1 = sinh(len);
      float c2 = cosh(len) - 1.0;
      return mat4(1.0) + c1 * m + c2 * m2;
    }
}

//-------------------------------------------------------

vec4 getRay(float fov, vec2 resolution, vec2 fragCoord){
  vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
  float z = 0.1;
  vec3 pPre = qtransform(cameraQuat, vec3(-xy,z));
  vec4 p = projectToHyperboloid(vec4(pPre, 1.0));
  return p;
}

//-------------------------------------------------------
//Raymarch functions
//-------------------------------------------------------
float unionSDF(float d1, float d2){
  return min(d1, d2);
}

float differenceSDF(float d1, float d2){
  return max(-d1, d2);
}

//Raymarch primitives

float sphereHSDF(vec4 samplePoint, vec4 center, float radius){
  return hypDistance(samplePoint, center) - radius;
}

float horosphereHSDF(vec4 samplePoint, vec4 lightPoint){
  return log(lorentzDot(samplePoint, lightPoint));
}

float geodesicPlaneHSDF(vec4 samplePoint, vec4 dualPoint, float offset){
  return asinh(lorentzDot(samplePoint, dualPoint)) - offset;
}

float geodesicCylinderHSDF(vec4 samplePoint, vec4 dualPoint1, vec4 dualPoint2, float radius){
  // defined by two perpendicular geodesic planes
  return asinh(sqrt(pow(lorentzDot(samplePoint, dualPoint1),2.0)+pow(lorentzDot(samplePoint, dualPoint2),2.0))) - radius;
}

float sceneHSDF(vec4 samplePoint){  /// for {4,3,6} edges
   samplePoint = abs(samplePoint);
   //now reflect until smallest xyz coord is z, so it will be close to the xy edge of cube
   if(samplePoint.x < samplePoint.z){
    samplePoint = vec4(samplePoint.z,samplePoint.y,samplePoint.x,samplePoint.w);
   }
   if(samplePoint.y < samplePoint.z){
    samplePoint = vec4(samplePoint.x,samplePoint.z,samplePoint.y,samplePoint.w);
   }

   // should precompute these orthonomal calculations
   vec4 dualPoint1 = lorentzNormalize(vec4(1.0/halfCubeWidthKlein,0.0,0.0,1.0));
   vec4 dualPoint2 = vec4(0.0,1.0/halfCubeWidthKlein,0.0,1.0);
   dualPoint2 = lorentzNormalize(dualPoint2 + lorentzDot(dualPoint2, dualPoint1) * dualPoint1); 
  
   float final = geodesicCylinderHSDF(samplePoint, dualPoint1, dualPoint2, 0.1);
   return final;
 }

// float sceneHSDF(vec4 samplePoint){  /// for {4,3,6}
//    float sphereInv = -sphereHSDF(samplePoint, ORIGIN, sphereRad);
//    float horosphere = horosphereHSDF(abs(samplePoint), horosphereSize*idealCubeCornerKlein);
//    float diff = differenceSDF(horosphere, sphereInv);
//    float final = differenceSDF(horosphere, sphereInv);
//    // float final = horosphere;
//    return final;
//  }

/*float sceneHSDF(vec4 samplePoint){   /// for {4,3,7}
  float sphereInv = -sphereHSDF(samplePoint, ORIGIN, sphereRad);
  vec4 dualPoint = projectToHyperboloid(vec4(halfCubeWidthKlein,halfCubeWidthKlein,halfCubeWidthKlein,1.0));
  float plane0 = geodesicPlaneHSDF(abs(samplePoint), dualPoint, planeOffset);
  float diff = differenceSDF(plane0, sphereInv);
  float final = diff;
  return final;
}*/

// float sceneHSDF(vec4 samplePoint){   /// draw sides of the cube fundamental domain
//   // float sphereInv = -sphereHSDF(samplePoint, ORIGIN, sphereRad);
//   // float horosphere = horosphereHSDF(abs(samplePoint), horosphereSize*idealCubeCornerKlein);
//   // float diff = differenceSDF(horosphere, sphereInv);
//   vec4 dualPoint0 = projectToHyperboloid(vec4(1.0/halfCubeWidthKlein,0.0,0.0,1.0));
//   vec4 dualPoint1 = projectToHyperboloid(vec4(0.0,1.0/halfCubeWidthKlein,0.0,1.0));
//   vec4 dualPoint2 = projectToHyperboloid(vec4(0.0,0.0,1.0/halfCubeWidthKlein,1.0));
//   float plane0 = geodesicPlaneHSDF(abs(samplePoint), dualPoint0, 0.0);
//   float plane1 = geodesicPlaneHSDF(abs(samplePoint), dualPoint1, 0.0);
//   float plane2 = geodesicPlaneHSDF(abs(samplePoint), dualPoint2, 0.0);
//   // float final = unionSDF(unionSDF(unionSDF(diff,plane0),plane1),plane2);
//   float final = unionSDF(unionSDF(plane0,plane1),plane2);
//   // float final = differenceSDF(horosphere, sphereInv);
//   return final;
// }


//-------------------------------------------------------


bool isOutsideCell(vec4 samplePoint, out mat4 fixMatrix){
  vec4 kleinSamplePoint = projectToKlein(samplePoint);
  if(kleinSamplePoint.x > halfCubeWidthKlein){
    fixMatrix = invGenerators[0];
    return true;
  }
  if(kleinSamplePoint.x < -halfCubeWidthKlein){
    fixMatrix = invGenerators[1];
    return true;
  }
  if(kleinSamplePoint.y > halfCubeWidthKlein){
    fixMatrix = invGenerators[2];
    return true;
  }
  if(kleinSamplePoint.y < -halfCubeWidthKlein){
    fixMatrix = invGenerators[3];
    return true;
  }
  if(kleinSamplePoint.z > halfCubeWidthKlein){
    fixMatrix = invGenerators[4];
    return true;
  }
  if(kleinSamplePoint.z < -halfCubeWidthKlein){
    fixMatrix = invGenerators[5];
    return true;
  }
  return false;
}

float raymarchDistance(vec4 rO, vec4 rD, float start, float end, out vec4 endPoint, out vec4 endRayTangentVector, out float tilingSteps){
  int fakeI = 0;
  float totalDepth = start;
  float localDepth = totalDepth;
  mat4 fixMatrix;
  for(int i = 0; i< MAX_MARCHING_STEPS; i++){
    if(fakeI >= maxSteps){
      //when we break its as if we reached our max marching steps
      break;
    }
    fakeI++;
    vec4 samplePoint = pointOnGeodesic(rO, rD, localDepth);
    if(isOutsideCell(samplePoint, fixMatrix)){
      tilingSteps++;
      vec4 newDirection = pointOnGeodesic(rO, rD, localDepth + 0.1); //forwards a bit
      rO = samplePoint*fixMatrix;
      newDirection *= fixMatrix;
      rO = projectToHyperboloid(rO);
      newDirection = projectToHyperboloid(newDirection);
      rD = vPrimeFromV(rO,newDirection);
      localDepth = start;
    }
    else{
      float dist = sceneHSDF(samplePoint);
      //float dist = unionSDF(sceneHSDF(samplePoint), sphereHSDF(samplePoint, rO*translateByVector(vec3(0,0,-0.2)), 0.1));
      if(dist < EPSILON){
        endPoint = samplePoint;
        endRayTangentVector = tangentVectorOnGeodesic(rO, rD, localDepth);
        return totalDepth;
      }
      totalDepth += dist;
      localDepth += dist;
      if(totalDepth >= end){
        endPoint = pointOnGeodesic(rO, rD, localDepth);
        endRayTangentVector = tangentVectorOnGeodesic(rO, rD, localDepth);
        return end;
      }
    }
  }
  endPoint = pointOnGeodesicAtInfinity(rO, rD);
  endRayTangentVector = tangentVectorOnGeodesic(rO, rD, localDepth);
  return end;
}

//COLORING FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++
vec4 estimateNormal(vec4 p, vec4 rO) { // normal vector is in tangent plane to hyperboloid at p
  // float denom = sqrt(1.0 + p.x*p.x + p.y*p.y + p.z*p.z);  // first, find basis for that tangent hyperplane
  float denom = p.w;
  vec4 basis_x = lorentzNormalize(vec4(denom,0.0,0.0,p.x));  // dw/dx = x/denom on hyperboloid
  vec4 basis_y = vec4(0.0,denom,0.0,p.y);  // dw/dy = y/denom
  vec4 basis_z = vec4(0.0,0.0,denom,p.z);  // dw/dz = z/denom  /// note that these are not orthonormal!
  basis_y = lorentzNormalize(basis_y + lorentzDot(basis_y, basis_x)*basis_x); // need to Gram Schmidt
  basis_z = lorentzNormalize(basis_z + lorentzDot(basis_z, basis_x)*basis_x + lorentzDot(basis_z, basis_y)*basis_y);
 return lorentzNormalize(
     basis_x * (sceneHSDF(projectToHyperboloid(p + EPSILON*basis_x)) - sceneHSDF(projectToHyperboloid(p - EPSILON*basis_x))) +
     basis_y * (sceneHSDF(projectToHyperboloid(p + EPSILON*basis_y)) - sceneHSDF(projectToHyperboloid(p - EPSILON*basis_y))) +
     basis_z * (sceneHSDF(projectToHyperboloid(p + EPSILON*basis_z)) - sceneHSDF(projectToHyperboloid(p - EPSILON*basis_z)))
 );
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void main(){
  vec4 endPoint = vec4(0.0,0.0,0.0,1.0);
  vec4 endRayTangentVector = vec4(0.0,0.0,0.0,0.0);
  float tilingSteps = 1.0;
  vec4 rayOrigin = vec4(0.0,0.0,0.0,1.0);
  //camera position must be translated in hyperboloid ------------------------
  //rayOrigin *= translateByVector(cameraPos);
  rayOrigin *= currentBoost;
  //generate direction then transform to hyperboloid ------------------------
  vec4 rayDirV = getRay(90.0, screenResolution, gl_FragCoord.xy);
  //rayDirV *= translateByVector(cameraPos);
  rayDirV *= currentBoost;
  vec4 rayDirVPrime = vPrimeFromV(rayOrigin, rayDirV);
  //get our raymarched distance back ------------------------
  float dist = raymarchDistance(rayOrigin, rayDirVPrime, MIN_DIST, MAX_DIST, endPoint, endRayTangentVector, tilingSteps);
  if((dist > MAX_DIST - EPSILON)||(tilingSteps >= float(MAX_MARCHING_STEPS) - 0.5)){
    //Didn't hit anything ------------------------
    vec4 pointAtInfinity = pointOnGeodesicAtInfinity(rayOrigin, rayDirVPrime);
    gl_FragColor = vec4(0.5*pointAtInfinity.xyz+vec3(0.5,0.5,0.5),1.0);
    return;
  }

  vec4 surfaceNormal = estimateNormal(endPoint, rayOrigin);
  float shineShade = lorentzDot(surfaceNormal, endRayTangentVector);
  float depthShade = max(1.0-dist/5.0, 0.0);
  float stepsShade = max(1.0-tilingSteps/3.0,0.0);
  // float comboShade = shineShade*depthShade;
  vec4 depthColor = vec4(depthShade,depthShade*0.65,0.1,1.0);
  // vec4 stepsColor = vec4(stepsShade,stepsShade,stepsShade,1.0);
  vec4 shineColor = vec4(shineShade,shineShade,shineShade,1.0);
  // vec4 comboColor = vec4(comboShade,comboShade,comboShade,1.0);
  // vec4 orange = vec4(1.0,0.65,0.1,1.0);
  // vec4 white = vec4(1.0,1.0,1.0,1.0);
  // vec4 normalColor = vec4(abs(normalize(projectToKlein(surfaceNormal).xyz)),1.0);
  //abs is needed to avoid inconsistencies in shading coming from different paths
  // vec4 endRayColor = vec4((normalize(projectToKlein(endRayTangentVector).xyz)),1.0);
  //to the same cube giving different orientations of the cube

  // gl_FragColor = 0.85 * depthColor + 0.15 * normalColor;
  // if(comboShade < 0.5){
  //   gl_FragColor = 2.0 * comboShade * orange;
  // }
  // else{
  //   gl_FragColor = 2.0*(comboShade-0.5)*white + (1.0 - 2.0*(comboShade-0.5))*orange;
  // }
  gl_FragColor = 0.5*depthColor + 0.5*shineColor;
  // gl_FragColor = shineColor;
  // gl_FragColor = 0.2*stepsColor + 0.8*normalColor;
  // gl_FragColor = normalColor;
  // gl_FragColor = endRayColor;
}
