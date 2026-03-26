import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

class NURBSLab {
    constructor() {
        this.container = document.getElementById('three-container');
        this.engine = new window.NURBS2D();
        
        this.setupScene();
        this.setupPatch();
        this.setupUI();
        this.render();

        window.addEventListener('resize', () => this.onWindowResize());
    }

    setupScene() {
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x020617);

        this.camera = new THREE.PerspectiveCamera(75, this.container.clientWidth / this.container.clientHeight, 0.1, 1000);
        this.camera.position.set(3, 3, 5);

        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;

        // Lights
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);

        const pointLight = new THREE.PointLight(0x3b82f6, 1, 100);
        pointLight.position.set(5, 5, 5);
        this.scene.add(pointLight);

        // Helpers
        const grid = new THREE.GridHelper(10, 10, 0x334155, 0x1e293b);
        grid.rotation.x = Math.PI / 2;
        this.scene.add(grid);
    }

    setupPatch() {
        // Initial Quadratic NURBS Surface Data
        this.patch = {
            p: 2,
            q: 2,
            U: [0, 0, 0, 1, 1, 1],
            V: [0, 0, 0, 1, 1, 1],
            weights: [
                [1, 1, 1],
                [1, 1, 1],
                [1, 1, 1]
            ],
            controlPoints: [
                [ {x: -1, y: -1, z: 0}, {x: 0, y: -1, z: 0}, {x: 1, y: -1, z: 0} ],
                [ {x: -1, y: 0, z: 0},  {x: 0, y: 0, z: 1.5},  {x: 1, y: 0, z: 0} ],
                [ {x: -1, y: 1, z: 0},  {x: 0, y: 1, z: 0}, {x: 1, y: 1, z: 0} ]
            ]
        };

        // Surface Mesh
        this.geometry = new THREE.PlaneGeometry(1, 1, 30, 30);
        this.material = new THREE.MeshPhongMaterial({
            color: 0x3b82f6,
            side: THREE.DoubleSide,
            wireframe: false,
            flatShading: false,
            shininess: 100
        });
        this.mesh = new THREE.Mesh(this.geometry, this.material);
        this.scene.add(this.mesh);

        // Control Points Visualization
        this.cpGroup = new THREE.Group();
        this.scene.add(this.cpGroup);
        
        this.updateSurface();
    }

    updateSurface() {
        const positions = this.geometry.attributes.position.array;
        const resolution = 30; // Matches PlaneGeometry segments
        
        for (let i = 0; i <= resolution; i++) {
            for (let j = 0; j <= resolution; j++) {
                const u = i / resolution;
                const v = j / resolution;
                
                const point = this.engine.evaluateSurface(this.patch, u, v);
                
                const index = (i * (resolution + 1) + j) * 3;
                positions[index] = point.x;
                positions[index + 1] = point.y;
                positions[index + 2] = point.z;
            }
        }
        
        this.geometry.attributes.position.needsUpdate = true;
        this.geometry.computeVertexNormals();

        // Update CP helpers
        this.cpGroup.clear();
        for (let i = 0; i < this.patch.controlPoints.length; i++) {
            for (let j = 0; j < this.patch.controlPoints[0].length; j++) {
                const cp = this.patch.controlPoints[i][j];
                const sphere = new THREE.Mesh(
                    new THREE.SphereGeometry(0.05),
                    new THREE.MeshBasicMaterial({ color: 0x10b981 })
                );
                sphere.position.set(cp.x, cp.y, cp.z);
                this.cpGroup.add(sphere);
            }
        }
    }

    setupUI() {
        const resU = document.getElementById('resU');
        const resV = document.getElementById('resV');
        const pHeight = document.getElementById('pHeight');
        const pTwist = document.getElementById('pTwist');

        resU.addEventListener('input', (e) => {
            const val = parseInt(e.target.value);
            document.getElementById('resU-val').innerText = val;
            this.rebuildSurface(val, parseInt(resV.value));
        });

        resV.addEventListener('input', (e) => {
            const val = parseInt(e.target.value);
            document.getElementById('resV-val').innerText = val;
            this.rebuildSurface(parseInt(resU.value), val);
        });

        pHeight.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('height-val').innerText = val.toFixed(1);
            this.patch.controlPoints[1][1].z = val; // Move center point
            this.updateSurface();
        });

        pTwist.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('twist-val').innerText = val.toFixed(1);
            this.patch.controlPoints[0][0].z = val;
            this.patch.controlPoints[2][2].z = -val;
            this.updateSurface();
        });
    }

    rebuildSurface(u, v) {
        this.scene.remove(this.mesh);
        this.geometry.dispose();
        this.geometry = new THREE.PlaneGeometry(1, 1, u, v);
        this.mesh = new THREE.Mesh(this.geometry, this.material);
        this.scene.add(this.mesh);
        this.updateSurface();
    }

    onWindowResize() {
        this.camera.aspect = this.container.clientWidth / this.container.clientHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
    }

    render() {
        requestAnimationFrame(() => this.render());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}

// Start lab
new NURBSLab();
