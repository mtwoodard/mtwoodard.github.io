BEGIN VERTEX
void main()
{
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position.xyz, 1.0);
}
END VERTEX

BEGIN FRAGMENT
const int MAX_MARCHING_STEPS = 127;

uniform int isStereo;
uniform vec2 screenResolution;
uniform mat4 currentBoost;
uniform float halfCubeWidthKlein;
uniform mat4 invGenerators[6];

//--------------------------------------------------------------------
// Hyperbolic Math Functions
//--------------------------------------------------------------------
float cosh(float x){
  float eX = exp(x);
  return (0.5 * (eX + 1.0/eX));
}

float acosh(float x){ //must be more than 1
  return log(x + sqrt(x*x-1.0));
}

float sinh(float x){
  float eX = exp(x);
  return (0.5 * (eX - 1.0/eX));
}

float asinh(float x){
  return log(x + sqrt(x*x+1.0));
}

//-------------------------------------------------------
// Hyperbolic Functions
//-------------------------------------------------------
float hypDot(vec4 u, vec4 v){
  return u.x*v.x + u.y*v.y + u.z*v.z - u.w*v.w; // hyp Dot
}
float hypNorm(vec4 v){
  return sqrt(abs(hypDot(v,v)));
}
vec4 hypNormalize(vec4 u, bool toTangent){
  return u/hypNorm(u);
}
float hypDistance(vec4 u, vec4 v){
  float bUV = -hypDot(u,v);
  return acosh(bUV);
}

//Given two positions find the unit tangent vector at the first that points to the second
vec4 hypDirection(vec4 u, vec4 v){
  vec4 w = v + hypDot(u,v)*u;
  return hypNormalize(w, true);
}

//calculate the new direction vector (v) for the continuation of the ray from the new ray origin (u)
//having moved by fix matrix
vec4 hypFixDirection(vec4 u, vec4 v, mat4 fixMatrix){
  return hypDirection(u, v*fixMatrix); 
}

vec4 projectToKlein(vec4 v){
  return v/v.w;
}

//-------------------------------------------------------
// Raymarch Function and Helpers
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

vec4 getRayPoint(vec2 resolution, vec2 fragCoord, bool isLeft){ //creates a point that our ray will go through
    if(isStereo == 1){
      resolution.x = resolution.x * 0.5;
      if(!isLeft) { fragCoord.x = fragCoord.x - resolution.x; }
    }
    vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
    float z = 0.1/tan(radians(90.0*0.5));
    vec4 p =  hypNormalize(vec4(xy,-z,1.0), false);
    return vec4(xy,-z,1.0);
}

void raymarch(vec4 rO, vec4 rD, out mat4 totalFixMatrix){
    //March Local Scene First
    for(int i = 0; i< MAX_MARCHING_STEPS; i++){

    }

    //March Global Scene
    for(int i = 0; i< MAX_MARCHING_STEPS; i++){


    }
}

void main(){
    vec4 rayOrigin = vec4(0,0,0,1);
        
    bool isLeft = gl_FragCoord.x/screenResolution.x <= 0.5;
    vec4 rayDirV = getRayPoint(screenResolution, gl_FragCoord.xy, isLeft);

    rayOrigin *= currentBoost;
    rayDirV *= currentBoost;

    vec4 rayDirVPrime = hypDirection(rayOrigin, rayDirV);

    //mat4 totalFixMatrix = mat4(1.0);
    //raymarch(rayOrigin, rayDirVPrime, totalFixMatrix);

    gl_FragColor = rayDirVPrime;
}

END FRAGMENT