/**
 * @author alteredq / http://alteredqualia.com/
 */

THREE.ShaderPass = function ( shader, textureID ) {

	THREE.Pass.call( this );

	this.textureID = ( textureID !== undefined ) ? textureID : "tDiffuse";

	if ( shader instanceof THREE.ShaderMaterial ) {

		this.uniforms = shader.uniforms;

		this.material = shader;

	} else if ( shader ) {

		this.uniforms = THREE.UniformsUtils.clone( shader.uniforms );

		this.material = new THREE.ShaderMaterial( {

			defines: Object.assign( {}, shader.defines ),
			uniforms: this.uniforms,
			vertexShader: shader.vertexShader,
			fragmentShader: shader.fragmentShader

		} );

	}

	this.camera = new THREE.OrthographicCamera( - 1, 1, 1, - 1, 0, 1 );
	this.scene = new THREE.Scene();

	this.quad = new THREE.Mesh( new THREE.PlaneBufferGeometry( 2, 2 ), null );
	this.quad.frustumCulled = false; // Avoid getting clipped
	this.scene.add( this.quad );

};

THREE.ShaderPass.prototype = Object.assign( Object.create( THREE.Pass.prototype ), {

	constructor: THREE.ShaderPass,

	render: function( renderer, writeBuffer, readBuffer, delta, maskActive ) {
		if(g_vr) 
			this.renderStereo(renderer, writeBuffer, readBuffer, delta, maskActive);
		else 
			this.renderMono(renderer, writeBuffer, readBuffer, delta, maskActive);
	},

	renderMono: function( renderer, writeBuffer, readBuffer, delta, maskActive ) {

		if ( this.uniforms[ this.textureID ] ) {
			this.uniforms[ this.textureID ].value = readBuffer.texture;
		}

		this.quad.material = this.material;

		if ( this.renderToScreen ) {
			renderer.render( this.scene, this.camera );
		} 
		else {
			renderer.render( this.scene, this.camera, writeBuffer, this.clear );
		}

	},

	renderStereo: function(renderer, writeBuffer, readBuffer, delta, maskActive){
		if ( this.uniforms[ this.textureID ] ) {
			this.uniforms[ this.textureID ].value = readBuffer.texture;
		}

		this.quad.material = this.material;

		var size = renderer.getSize();
		var rendererWidth = size.width;
		var rendererHeight = size.height;
		var eyeDivisionLine = rendererWidth/2;

		renderer.setScissorTest(true);

		if(this.renderToScreen){
			//render left 
			if(this.uniforms.isStereo)
				this.uniforms.isStereo.value = -1;
			renderer.setViewport(0, 0, eyeDivisionLine, rendererHeight);
			renderer.setScissor(0, 0, eyeDivisionLine, rendererHeight);
			renderer.render(this.scene, this.camera);

			//render right eye
			if(this.uniforms.isStereo)
				this.uniforms.isStereo.value = 1;
			renderer.setViewport(eyeDivisionLine, 0, eyeDivisionLine, rendererHeight);
			renderer.setScissor(eyeDivisionLine, 0, eyeDivisionLine, rendererHeight);
			renderer.render(this.scene, this.camera);
		}
		else{
			//render left eye
			if(this.uniforms.isStereo)
				this.uniforms.isStereo.value = -1;
			renderer.setViewport(0, 0, eyeDivisionLine, rendererHeight);
			renderer.setScissor(0, 0, eyeDivisionLine, rendererHeight);
			renderer.render(this.scene, this.camera, writeBuffer, this.clear );

			//render right eye
			if(this.uniforms.isStereo)
				this.uniforms.isStereo.value = 1;
			renderer.setViewport(eyeDivisionLine, 0, eyeDivisionLine, rendererHeight);
			renderer.setScissor(eyeDivisionLine, 0, eyeDivisionLine, rendererHeight);
			renderer.render(this.scene, this.camera, writeBuffer, this.clear );
		}

	}

} );
