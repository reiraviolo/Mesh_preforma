"""Malla de cuadriláteros de un tee con gmsh (OpenCASCADE)."""
import numpy as np, gmsh, sys

def malla_tee(OD_1=200, OD_2=200, H_1=200, H_2=200, H_3=200, Alpha_2=90, tam=6.28, salida=None):
    R1, R2 = OD_1/2, OD_2/2
    a = np.radians(Alpha_2); eje = np.array([np.cos(a), 0.0, np.sin(a)])
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("tee")
    o = gmsh.model.occ
    directo = o.addCylinder(-H_1, 0, 0, H_1 + H_3, 0, 0, R1)
    deriv = o.addCylinder(0, 0, 0, *(eje * H_2), R2)
    o.fuse([(3, directo)], [(3, deriv)])
    o.synchronize()
    # se descartan las 3 tapas planas: la superficie de la cáscara es el resto
    caras = []
    for dim, tag in gmsh.model.getEntities(2):
        tipo = gmsh.model.getType(dim, tag)
        if "Plane" not in tipo:
            caras.append(tag)
    gmsh.model.addPhysicalGroup(2, caras, name="cascara")
    gmsh.option.setNumber("Mesh.MeshSizeMin", tam * 0.95)
    gmsh.option.setNumber("Mesh.MeshSizeMax", tam * 1.05)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)          # Frontal-Delaunay para cuadriláteros
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # blossom (no soporta full-quad en superficies periódicas)
    gmsh.option.setNumber("Mesh.Smoothing", 20)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")
    etiquetas, coords, _ = gmsh.model.mesh.getNodes()
    orden = np.argsort(etiquetas)
    nodos = coords.reshape(-1, 3)[orden]
    mapa = np.zeros(int(etiquetas.max()) + 1, dtype=int); mapa[etiquetas[orden]] = np.arange(len(etiquetas))
    quads, tris = [], []
    for dim, tag in gmsh.model.getEntities(2):
        tipos, _, nodos_el = gmsh.model.mesh.getElements(2, tag)
        for t, nn in zip(tipos, nodos_el):
            if t == 3: quads.append(mapa[nn.reshape(-1, 4)])
            elif t == 2: tris.append(mapa[nn.reshape(-1, 3)])
    quads = np.vstack(quads) if quads else np.zeros((0, 4), int)
    tris = np.vstack(tris) if tris else np.zeros((0, 3), int)
    if salida: gmsh.write(salida)
    gmsh.finalize()
    return nodos, quads, tris

if __name__ == "__main__":
    for al in [30, 60, 90]:
        n, q, t = malla_tee(Alpha_2=al, H_2=500, tam=6.28)
        p = n[q]
        h = np.sqrt(0.5*np.linalg.norm(np.cross(p[:,2]-p[:,0], p[:,3]-p[:,1]), axis=1))
        forma = 4*0.5*np.linalg.norm(np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1]),axis=1)/np.sum(
            [np.linalg.norm(p[:,(i+1)%4]-p[:,i],axis=1)**2 for i in range(4)],axis=0)
        nr = np.cross(p[:,2]-p[:,0], p[:,3]-p[:,1]); nr /= np.linalg.norm(nr,axis=1,keepdims=True)
        sg = np.array([np.einsum('ij,ij->i', np.cross(p[:,(i+1)%4]-p[:,i], p[:,(i-1)%4]-p[:,i]), nr) for i in range(4)])
        print("%3d°: %5d cuads + %d triangulos | forma med %.3f min %.3f | h p5 %.2f med %.2f p95 %.2f | disp %.2f | plegados %d"
              % (al, len(q), len(t), forma.mean(), forma.min(), np.percentile(h,5), np.median(h), np.percentile(h,95),
                 np.percentile(h,95)/np.percentile(h,5), int((sg.min(0)<0).sum())))
